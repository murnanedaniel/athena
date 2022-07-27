# Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration

__doc__ = """This service runs in the online Athena event display threads. It
manages the distribution of incoming events to the right event display streams.
In order to do that it connects to the different Atlantis and VP1 event output
algorithms/tools and manipulates their job options during the run.
"""
__version__ = "0.1"
__author__  = "Eric Jansen <eric.jansen@cern.ch>"

import os, grp, stat, random
from ipc import IPCPartition
from ispy import ISInfoDictionary, ISInfoAny, ISObject
from AthenaPython import PyAthena
from AthenaPython.PyAthena import StatusCode
from PyAnalysisCore import PyEventTools
from EventDisplaysOnline import EventUtils

class OnlineEventDisplaysSvc( PyAthena.Svc ):

	def __init__(self, name = 'OnlineEventDisplaysSvc', **kw):
		kw['name'] = name
		self.output = kw.get('OutputDirectory')
		self.maxevents = kw.get('MaxEvents')
		self.public = kw.get('Public')
		self.projecttags = kw.get('ProjectTags')

		self.zpgid = None
		self.partition = None

		self.StreamToFileTool = None
		self.StreamToServerTool = None
		self.VP1EventProducer = None

		self.run = 0
		self.event = 0
		self.stream = ''
		self.directory = ''

		super(OnlineEventDisplaysSvc,self).__init__(**kw)

	def initialize(self):
		incidentSvc = PyAthena.py_svc('IncidentSvc', iface='IIncidentSvc')
		incidentSvc.addListener(self, 'BeginEvent')
		incidentSvc.addListener(self, 'StoreCleared')

		self.partition = IPCPartition('ATLAS')
		self.dict = ISInfoDictionary(self.partition)
		self.zpgid = grp.getgrnam("zp").gr_gid
		return StatusCode.Success

	def finalize(self):
		return StatusCode.Success

	def getJobOptions(self):
		self.StreamToFileTool = PyAthena.py_tool('JiveXML::StreamToFileTool', iface='IProperty')
		self.StreamToServerTool = PyAthena.py_tool('JiveXML::StreamToServerTool', iface='IProperty')

		# Some casting magic is needed to access algorithm properties
		from GaudiPython.Bindings import gbl, InterfaceCast
		vp1alg = PyAthena.py_alg('VP1EventProd')
		self.VP1EventProducer = InterfaceCast(gbl.IProperty).cast(vp1alg)

		self.msg.info("StreamToFileTool: %s", self.StreamToFileTool)
		self.msg.info("StreamToServerTool: %s", self.StreamToServerTool)
		self.msg.info("VP1EventProducer: %s", self.VP1EventProducer)

	def beginEvent(self):
		if not (self.StreamToFileTool and self.StreamToServerTool and self.VP1EventProducer):
			self.getJobOptions()

		try:
			eventInfo = PyEventTools.getEventInfo('EventInfo')
		except LookupError as err:
			self.msg.error("Could not retrieve EventInfo: %s", err)
			return StatusCode.Recoverable

		try:
			# Read event info
			self.run = eventInfo.runNumber()
			self.event = eventInfo.eventNumber()

			# Retrieve trigger info
			streamTags = eventInfo.streamTags()
		except Exception as err:
			self.msg.error("Exception occured while reading event/trigger info: %s", err)
			return StatusCode.Recoverable

		# Retrieve the physics stream names from the trigger info
		streams = []
		for tag in streamTags:
			
			## Added 08/03/15 - sjiggins - Wanted to see if events had physics stream tag or was from random trigger
			if tag.name():
				self.msg.debug("Event %d/%d has the corresponding streamTags: %s", self.run, self.event, tag.type())
			##################################################################################
			if tag.type() == 'express' and tag.name():
				streams += [tag.type()+'_'+tag.name()]

			if tag.type() == 'physics' and tag.name():
				streams += [tag.type()+'_'+tag.name()]

		# Add special streams to the list Public
		try:
			for stream in streams:
				if stream in self.public:
					ready4physics = ISInfoAny()
					self.dict.getValue('RunParams.Ready4Physics', ready4physics)
					print("Ready for physics: %s " % ready4physics.get())
					runparams = ISObject(self.partition, 'RunParams.RunParams','RunParams')
					runparams.checkout()
					physicsReady = ISObject(self.partition, 'RunParams.Ready4Physics','Ready4PhysicsInfo')
					physicsReady.checkout()
					print("Ready for physics: %r" % (physicsReady.ready4physics))
					#if ready4physics.get() and physicsReady.ready4physics and runparams.T0_project_tag in self.projecttags:
					if physicsReady.ready4physics and runparams.T0_project_tag in self.projecttags:
						streams += ['Public']
					else:
						self.msg.debug("RunParams.Ready4Physics is not set, run number is not set, or T0_project_tag is not set to any of %s", ", ".join(self.projecttags))
						break
		except Exception as err:
			self.msg.error("Exception occured while reading RunParams.Ready4Physics: %s", err)

		# Randomize list of streams
		random.shuffle(streams)
		self.msg.debug("Event %d/%d has event display stream tags: %s", self.run, self.event, ", ".join(streams))

		## for beam splash, give priority to MinBias
		#if 'physics_MinBias' in streams:
		#	streams = ['physics_MinBias']
		#	self.msg.debug("Modified stream tag: %s", ", ".join(streams))

		# Start from the beginning and send the event to the first stream that passes our directory checks
		self.directory = ''
		for self.stream in streams:
			self.directory = "%s/%s" % (self.output, self.stream)
			if os.access(self.directory, os.F_OK):
				if os.path.isdir(self.directory) and os.access(self.directory, os.W_OK):
					self.msg.debug("Going to write file to existing directory: %s", self.directory)
					if os.stat(self.directory).st_gid != self.zpgid:
						self.msg.debug("Setting group to 'zp' for directory: %s", self.directory)
						os.chown(self.directory, -1, self.zpgid)
					break
				else:
					self.msg.warning("Directory \'%s\' is not usable, trying next alternative", self.directory)
					self.directory = ''
			else:
				try:
					os.mkdir(self.directory)
					os.chmod(self.directory, stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)
					os.chown(self.directory, -1, self.zpgid)
					self.msg.info("Created output directory \'%s\' for stream \'%s\'", self.directory, self.stream)
					break
				except OSError as err:
					self.msg.warning("Failed to create output directory \'%s\' for stream \'%s\': %s", (self.directory, self.stream, err.strerror))
					self.directory = ''

		# Check if a suitable directory was found
		if self.directory:
			self.msg.debug("Event %d/%d will be streamed to: %s", self.run, self.event, self.stream)
		else:
			# This event is hopelessly lost, send StatusCode.Recoverable in an attempt to abort.
			# But if Athena chooses to ignore that, set the output to the "Unknown" trashcan stream.
			self.stream = 'Unknown'
			self.directory = "%s/.Unknown" % self.output
		try:
			# Set output stream for JiveXML event streaming (file and server)
			self.StreamToServerTool.getProperty('StreamName').setValue("%s" % self.stream)
			self.StreamToFileTool.getProperty('FileNamePrefix').setValue("%s/JiveXML" % self.directory)

		except Exception as err:
			self.msg.error("Exception occured while setting job options: %s", err)
			return StatusCode.Failure

		if not self.directory:
			return StatusCode.Recoverable

	def endEvent(self):
		# VP1 copied its ESD file of the previous event to the destination directory.
		# So we set the directory of the current event in endEvent to pass it to the next event.
		if self.VP1EventProducer:
			self.VP1EventProducer.getProperty('DestinationDirectory').setValue(self.directory)

		# Prune events and make index file for atlas-live.cern.ch
		if self.directory:
			# If VP1 event producer is missing, skip the pair check to make cleanup/sync work
			if self.VP1EventProducer:
				EventUtils.cleanDirectory(self.msg, self.directory, self.maxevents, checkpair=True)
			else:
				EventUtils.cleanDirectory(self.msg, self.directory, self.maxevents, checkpair=False)

			## for beam splashes, disable prepareFilesForTransfer() in EventUtils.py
			## Activate these lines to zip and transfer every JiveXML file.
			#event = "%05d" % int(self.event)
			#EventUtils.zipXMLFile(self.msg, self.directory, f'JiveXML_{self.run}_{event}.xml')

		# And cleanup the variables
		self.run = 0
		self.event = 0
		self.stream = ''
		self.directory = ''

	def handle(self, incident):
		self.msg.verbose("Received incident %s from %s" % (incident.type(), incident.source()))

		# Event and trigger info present, decide stream
		if incident.type() == 'BeginEvent' and incident.source() == 'BeginIncFiringAlg':
			self.beginEvent()

		# VP1 writes its file at EndEvent, so we can do cleanup at StoreCleared
		if incident.type() == 'StoreCleared' and incident.source() == 'StoreGateSvc':
			self.endEvent()

