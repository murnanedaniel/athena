#ifndef QATDATAANALYSIS_OPTPARSE_H
#define QATDATAANALYSIS_OPTPARSE_H
class HistogramManager;
#include <string>
#include <vector>


// 
struct HIOZeroToOne {

  //---------------Simple Struct, few public variables--------------//
  HistogramManager *output;                                         //
  bool              verbose;                                        //
  void optParse(int & argc, char ** & argv);                        //
  //---------------End importatant part-----------------------------//

  HIOZeroToOne(const std::string & driver="");
  ~HIOZeroToOne();


  private:

  HIOZeroToOne & operator =( const HIOZeroToOne &);
  HIOZeroToOne(const HIOZeroToOne &);
  class Clockwork;
  Clockwork *m_c;

};

// 
struct HIOOneToOne {

  //---------------Simple Struct, few public variables--------------//
  const HistogramManager *input;                                    //
  HistogramManager *output;                                         //
  bool              verbose;                                        //
  void optParse(int & argc, char ** & argv);                        //
  //---------------End importatant part-----------------------------//

  HIOOneToOne(const std::string & driver="");
  ~HIOOneToOne();


  private:

  HIOOneToOne & operator =( const HIOOneToOne &);
  HIOOneToOne(const HIOOneToOne &);
  class Clockwork;
  Clockwork *m_c;

};

// 
struct HIOOneToZero {

  //---------------Simple Struct, few public variables--------------//
  const HistogramManager *input;                                    //
  bool              verbose;                                        //
  void optParse(int & argc, char ** & argv);                        //
  //---------------End importatant part-----------------------------//

  HIOOneToZero(const std::string & driver="");
  ~HIOOneToZero();


  private:

  HIOOneToZero & operator =( const HIOOneToZero &);
  HIOOneToZero(const HIOOneToZero &);
  class Clockwork;
  Clockwork *m_c;

};

// 
struct HIONToOne {

  //---------------Simple Struct, few public variables--------------//
  std::vector<const HistogramManager *> input;                      //
  HistogramManager *output;                                         //
  bool              verbose;                                        //
  void optParse(int & argc, char ** & argv);                        //
  //---------------End importatant part-----------------------------//

  HIONToOne(const std::string & driver="");
  ~HIONToOne();


  private:

  HIONToOne & operator =( const HIONToOne &);
  HIONToOne(const HIONToOne &);

  class Clockwork;
  Clockwork *m_c;

};
// 
struct HIONToZero {

  //---------------Simple Struct, few public variables--------------//
  std::vector<const HistogramManager *> input;                      //
  bool              verbose;                                        //
  void optParse(int & argc, char ** & argv);                        //
  //---------------End importatant part-----------------------------//

  HIONToZero(const std::string & driver="");
  ~HIONToZero();


  private:

  HIONToZero & operator =( const HIONToZero &);
  HIONToZero(const HIONToZero &);
  class Clockwork;
  Clockwork *m_c;

};


class NumericInput {
  
 public:

  NumericInput();
  ~NumericInput();

  // Declare a parameter (name, doc, value)
  void declare(const std::string & name, const std::string & doc, double val);

  // Then override with user input:
  void optParse(int & argc, char ** & argv);

  // Print the list of parameters:
  std::string usage() const;

  // Get the value of a parameter, by name:
  double getByName  (const std::string & name) const;

 private:
 
  NumericInput & operator =( const NumericInput &);
  NumericInput(const NumericInput &);

  class Clockwork;
  Clockwork *m_c;

};



#endif
