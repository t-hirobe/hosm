/*! @brief
 * @date 2010/05/06  1:57:18
 * $Lastupdate: 2010/05/07  3:36:48 $
 */

#include "All.h"
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>

using namespace std;

namespace flow {
namespace exception {

void All::addInfo(const string& fileInfo, const string& msg)
{
  string dat =
    (boost::format("%s\n\t%s") % fileInfo % msg).str();
  messages_.push_back(dat);
}

string All::info() const
{
  ostringstream stream;
  BOOST_FOREACH( string msg, messages_) {
    stream << "#" <<  msg << endl;
  }
  return stream.str();
}



} // exception
} // flow

