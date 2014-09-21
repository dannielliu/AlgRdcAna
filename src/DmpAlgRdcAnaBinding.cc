#include <boost/python.hpp>
#include "DmpAlgRdcAna.h"

BOOST_PYTHON_MODULE(libDmpAlgRdcAna){
  using namespace boost::python;

  class_<DmpAlgRdcAna,boost::noncopyable,bases<DmpVAlg> >("DmpAlgRdcAna",init<>());
}
