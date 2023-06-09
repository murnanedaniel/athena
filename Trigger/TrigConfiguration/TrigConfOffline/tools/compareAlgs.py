#!/usr/bin/env python
###########################################################################
#    Copyright (C) 2009 by Miroslav Nozicka                                      
#    <nozicka@mail.desy.de>                                                             
#
# Copyright: See COPYING file that comes with this distribution
#
###########################################################################

import getopt, sys

def usage():
  """Prints a help message"""
  print "Usage: %s [database-options]" % \
        sys.argv[0].split('/')[-1]
  print "  -h|-?|--help                Issues this help message"
  print " Database options may be:"
  print "  -p|--password               The user password to the DB"
  print "  -u|--user <string>          The name of the user in the DB"
  print "  -n|--name <string>          The name of the DB inside the server"
  print "  -m|--host <string>          The name of the host where the DB server is running"
  print "  -sm1|--sm1 <string>         Super Master Key of the first setup"
  print "  -sm2|--sm2 <string>         Super Master Key of the second setup"


if __name__ == '__main__':
  import cx_Oracle
  print dir(cx_Oracle)
  print cx_Oracle.apilevel
  print cx_Oracle.version
  print cx_Oracle.SYSDBA
  
  db_user='atlas_trig_nozicka'
  db_passwd='*******'
  db_host='devdb10'
  db_name = ''
  sm1=-1
  sm2=-1

  short_opt = "h?p:u:n:m:sm1:sm2"
  long_opt = ['help', 'password=', 'user=', 'name=', 'host=','sm1=','sm2=']
  
  if len(sys.argv) == 1:
      print "Missing arguments"
      usage()
      sys.exit(1)
      
  #these are common bootstrap procedures for option processing
  try:
      opts, args = getopt.getopt(sys.argv[1:], short_opt, long_opt)
  except getopt.GetoptError, exc:
      print '\nError: '+exc.msg+'\n'
      usage()
      sys.exit(1)

  # The defaults
    
  #Read the options
  for o, a in opts:
    if o in ['-h', '-?', '--help']:
      usage()
      sys.exit(0)
    if o in ['-p', '--password']: db_passwd = a
    if o in ['-u', '--user']:     db_user = a
    if o in ['-n', '--name']:     db_name = a
    if o in ['-m', '--host']:     db_host = a
    if o in ['-sm1', '--sm1']:    sm1 = a
    if o in ['-sm2', '--sm2']:    sm2 = a
  
  connection = cx_Oracle.Connection(db_user, db_passwd, db_host)
  cursor = connection.cursor()
  print cursor
  
  if db_name=='' : db_name=db_user
  prepend =''
  if db_name != db_user and db_name != '':
    prepend = '%s.' % (db_name)
  
  table_names = []
  query =  " SELECT %sall_objects.object_name " % (prepend)
  query += " FROM %sall_objects " % (prepend)
  query += " WHERE %sall_objects.object_type='TABLE' " % (prepend)
  query += " AND %sall_objects.owner='%s' " % (prepend,db_name.upper())
  cursor.execute(query)
  result = cursor.fetchall()
  for column in result :
    table_names.append(column[0].upper())
  
  print len(table_names),'tables available'
  trigger_keys= []
  # Get super Master keys
  if 'SUPER_MASTER_TABLE' not in table_names :
    print 'SUPER_MASTER_TABLE not found in database'
    cursor.close()
    sys.exit(1)
  else :
    query =  " SELECT %sSUPER_MASTER_TABLE.SMT_ID "  % (prepend)
    query += " FROM %sSUPER_MASTER_TABLE " % (prepend)
    cursor.execute(query)
    result = cursor.fetchall()
    for column in result :
      trigger_keys.append([column[0],[],[]])
  
  if 'HLT_PRESCALE_SET' not in table_names :
    print 'HLT_PRESCALE_SET not found in database'
  else :
    for i in range(len(trigger_keys)) :
      query  = " SELECT %sHLT_PRESCALE_SET.HPS_ID " % (prepend)
      query += " FROM %sSUPER_MASTER_TABLE " % (prepend)
      query += " JOIN %sHLT_MASTER_TABLE ON (%sSUPER_MASTER_TABLE.SMT_HLT_MASTER_TABLE_ID = %sHLT_MASTER_TABLE.HMT_ID) " % (prepend, prepend, prepend)
      query += " JOIN %sHLT_TM_TO_PS ON (%sHLT_MASTER_TABLE.HMT_TRIGGER_MENU_ID = %sHLT_TM_TO_PS.HTM2PS_TRIGGER_MENU_ID) " % (prepend, prepend, prepend)
      query += " JOIN %sHLT_PRESCALE_SET ON (%sHLT_PRESCALE_SET.HPS_ID = %sHLT_TM_TO_PS.HTM2PS_PRESCALE_SET_ID) " % (prepend, prepend, prepend)
      query += " WHERE %sSUPER_MASTER_TABLE.SMT_ID=:smt_id " % (prepend)
      cursor.execute(query, smt_id=trigger_keys[i][0])
      result = cursor.fetchall()
      for column in result :
        trigger_keys[i][1].append(column[0])
        
  if 'L1_PRESCALE_SET' not in table_names :
    print 'L1_PRESCALE_SET not found in database'
  else :
    for i in range(len(trigger_keys)) :
      query  = " SELECT %sL1_PRESCALE_SET.L1PS_ID " % (prepend)
      query += " FROM %sSUPER_MASTER_TABLE " % (prepend)
      query += " JOIN %sL1_MASTER_TABLE ON (%sSUPER_MASTER_TABLE.SMT_L1_MASTER_TABLE_ID = %sL1_MASTER_TABLE.L1MT_ID) " % (prepend, prepend, prepend)
      query += " JOIN %sL1_TM_TO_PS ON (%sL1_MASTER_TABLE.L1MT_TRIGGER_MENU_ID = %sL1_TM_TO_PS.L1TM2PS_TRIGGER_MENU_ID) " % (prepend, prepend, prepend)
      query += " JOIN %sL1_PRESCALE_SET ON (%sL1_PRESCALE_SET.L1PS_ID = %sL1_TM_TO_PS.L1TM2PS_PRESCALE_SET_ID) " % (prepend, prepend, prepend)
      query += " WHERE %sSUPER_MASTER_TABLE.SMT_ID=:smt_id " % (prepend)
      cursor.execute(query, smt_id=trigger_keys[i][0])
      result = cursor.fetchall()
      for column in result :
        trigger_keys[i][2].append(column[0])
  
  if len(trigger_keys) > 0 :
    print 'Available trigger keys:'
    print 'SMKey','\t','HLTPSK','\t','LVL1PSK'
  for keys in trigger_keys :
    print keys[0],'\t',keys[1],'\t',keys[2]
  
  cursor.close()
  
  
  
