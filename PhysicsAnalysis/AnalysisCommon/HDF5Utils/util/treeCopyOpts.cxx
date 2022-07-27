/*
Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "treeCopyOpts.h"
#include <iostream>

#include <boost/program_options.hpp>

namespace H5Utils {

  AppOpts getTreeCopyOpts(int argc, char* argv[])
  {
    namespace po = boost::program_options;
    AppOpts app;
    std::string usage = "usage: " + std::string(argv[0]) + " <files>..."
      + " -o <output> [-h] [opts...]\n";
    po::options_description opt(usage + "\nConvert a root tree to HDF5");
    opt.add_options()
      ("in-file",
       po::value(&app.file.in)->required()->multitoken(),
       "input file name")
      ("out-file,o",
       po::value(&app.file.out)->required(),
       "output file name")
      ("tree-name,t",
       po::value(&app.file.tree)->default_value("", "found"),
       "tree to use, use whatever is there by default (or crash if multiple)")
      ("help,h", "Print help messages")
      ("branch-regex,r",
       po::value(&app.tree.branch_regex)->default_value(""),
       "regex to filter branches")
      ("vector-lengths,l",
       po::value(&app.tree.vector_lengths)->multitoken()->value_name("args..."),
       "max size of vectors to write")
      ("verbose,v",
       po::bool_switch(&app.tree.verbose),
       "print branches copied")
      ("n-entries,n",
       po::value(&app.tree.n_entries)->default_value(0, "all")->implicit_value(1),
       "number of entries to copy")
      ("chunk-size,c",
       po::value(&app.tree.chunk_size)->default_value(CHUNK_SIZE),
       "chunk size in HDF5 file")
      ("selection,s",
       po::value(&app.tree.selection)->default_value(""),
       "selection string applied to ntuples")
      ("print-interval,p",
       po::value(&app.tree.print_interval)->default_value(0, "never")->implicit_value(-1, "1%"),
       "print progress")

      ;
    po::positional_options_description pos_opts;
    pos_opts.add("in-file", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(opt)
                .positional(pos_opts).run(), vm);
      if ( vm.count("help") ) {
        std::cout << opt << std::endl;
        app.exit_code = 1;
      }
      po::notify(vm);
    } catch (po::error& err) {
      std::cerr << usage << "ERROR: " << err.what() << std::endl;
      app.exit_code = 1;
    }
    return app;
  }

}
