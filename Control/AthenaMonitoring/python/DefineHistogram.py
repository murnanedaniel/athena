def defineHistogram(varname, type='TH1F', path='EXPERT',
                    title='Unspecified_title_for_the_histogram_is_truly_annoying,_because_this_default_is_long;unspecified_label;unspecified_label',
                    xbins=100, xmin=0, xmax=1,
                    ybins=None, ymin=None, ymax=None, zmin=None, zmax=None, opt='', labels=None):
    """ Generates the histogram definition string which is digestable by TrigGenericMonitoringTool.


    Note that defaults are compleetly unreasonable. The tile is intentionally made annoying.
    For histogras types supported and allowed paths have a look at TrigGenericMonitoringTool documentation.
    """
    coded = "%s, %s, %s, %s, %d, %f, %f" % (path, type, varname, title, xbins, xmin, xmax)
    if ybins is not None:
        coded += ",%d, %f, %f" % (ybins, ymin, ymax)
        if zmin is not None:
            coded += ", %f, %f" % (zmin, zmax)
    if ybins is None and ymin is not None:
        coded += ", %f, %f" % (ymin, ymax)
    if labels is not None:
        labels = labels.strip()   # remove spurious white-spaces
        if len(labels)>0:
            if labels[-1]!=':': labels += ':'  # C++ parser expects at least one ":"
            coded += ",%s " % labels
    coded += ", %s" % opt

    return coded
