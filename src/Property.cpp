#include "Property.h"
#include "exception/All.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/exception/all.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <typeinfo>

namespace flow {
    namespace util {

        using namespace std;

        typedef boost::numeric::ublas::vector<double> UbVectorD;
/// Default constructor
        Property::Property() {}

        Property::Property(const string& fileName)
            : fileName_(fileName)
        {
            try {
                read(fileName_);
            } catch( flow::exception::IO& e) {
                string msg = "";
                flow::exception::IO e;
                MACRO_THROW_EXCEPTION( e, msg);
            }
        }

/// Default destructor
        Property::~Property() {}

        void Property::read(const string& fileName)
        {
            fileName_ = fileName;
            ifstream in;
            // �㳰�����
            in.exceptions(fstream::badbit | fstream::failbit);
            try {
                in.open(fileName.c_str());
                setData_(in);
            } catch(ios_base::failure& e0) {
                string msg = (boost::format("file I/O error %s") % fileName).str();
                flow::exception::IO e;
                MACRO_THROW_EXCEPTION( e, msg) e;
            }
        }

        void Property::setData_(ifstream& in)
        {
            string line;
            vector<string> items;
            try {
                while( getline(in, line)) {
                    // ��Ƭ��#�ʤ饹���å�
                    if( boost::starts_with(line, "#"))
                        continue;
                    boost::algorithm::split(items, line, boost::is_any_of("="));
                    // 1�Ԥ�ɬ��"="�ϰ��
                    if(items.size() != 2)
                        continue;
                    dataMap_.insert
                        ( make_pair
                          (
                              boost::trim_copy(items[0]),
                              boost::trim_copy(items[1])
                              )
                            );
                }
            } catch (ios_base::failure& e0) {
                if(!in.eof()) {
                    string msg = "failed to read file";
                    flow::exception::IO e;
                    MACRO_THROW_EXCEPTION( e, msg) e;
                }
            }
        }

        template<typename Type>
        Type Property::get(const string& key) const
        {
            if(!exist_(key)) {
                string msg =
                    (boost::format("key \"%s\" not found in property file: %s")
                     % key % fileName_).str();
                flow::exception::Parse e;
                MACRO_THROW_EXCEPTION(e, msg) e;
            }

            Type val;
            try {
                val = cast_<Type>( dataMap_.find(key)->second );
            } catch (flow::exception::Parse& e) {
                string msg = (boost::format("key : \"%s\"") % key).str();
                MACRO_THROW_EXCEPTION(e, msg);
            }
            return val;
        }
/*
        template<>
        UbVectorD Property::get(const string& key, const UbVectorD& defaultVal) const
        {
            if(!exist_(key))
                return defaultVal;

            UbVectorD val;
            try {
                val = cast_<UbVectorD>( dataMap_.find(key)->second );
            } catch (flow::exception::Parse& e) {
                string msg = (boost::format("key : \"%s\"") % key).str();
                MACRO_THROW_EXCEPTION(e, msg);
            }
            return val;
        }
*/
        template<typename Type>
        Type Property::get(const string& key, const Type& defaultVal) const
        {
            if(!exist_(key))
                return defaultVal;

            Type val;
            try {
                val = cast_<Type>( dataMap_.find(key)->second );
            } catch (flow::exception::Parse& e) {
                string msg = (boost::format("key : \"%s\"") % key).str();
                MACRO_THROW_EXCEPTION(e, msg);
            }
            return val;
        }

        void Property::output() const
        {
            pair<string, string> item;
            BOOST_FOREACH(item, dataMap_) {
                cerr << item.first << " = " << item.second << endl;
            }
        }

        template<typename Type>
        Type Property::cast_(const string& valStr) const
        {
            Type val;
            try {
                val = boost::lexical_cast<Type>( valStr );
            } catch(boost::bad_lexical_cast& e0) {
                string msg =
                    (boost::format("value \"%s\" is inconpatible with type of %s")
                     % valStr % string(typeid(Type).name())).str();
                flow::exception::Parse e;
                MACRO_THROW_EXCEPTION(e, msg) e;
            }
            return val;
        }

        template<>
        bool Property::cast_<bool>(const string& valStr) const
        {
            string val = boost::to_lower_copy(valStr);
            if(val == "true")
                return true;
            if(val == "false")
                return false;

            // true(TRUE) false(FALSE) �ʳ��Ǥ��㳰���ꤲ��
            string msg = "true(TREU) or false(FALSE) is allowed as bool type value";
            flow::exception::Parse e;
            MACRO_THROW_EXCEPTION(e, msg) e;
        }

        template<>
        UbVectorD Property::cast_<UbVectorD>(const string& valStr) const
        {
            string val = boost::to_lower_copy(valStr);
// valStr ��space��split
            UbVectorD vec;
            //vec�����ǿ���split��ʬ�䤷�����ˤ��碌��
            try {
                // ���Ǥ򥻥å�

                return vec;
            } catch(boost::bad_lexical_cast& e0) {
                flow::exception::Parse e;
                MACRO_THROW_EXCEPTION(e, "") e;
            }
        }

        template<>
        UbVectorD Property::get(const string& key, const UbVectorD& defaultVal) const
        {
            if(!exist_(key))
                return defaultVal;

            UbVectorD val;
            try {
                val = cast_<UbVectorD>( dataMap_.find(key)->second );
            } catch (flow::exception::Parse& e) {
                string msg = (boost::format("key : \"%s\"") % key).str();
                MACRO_THROW_EXCEPTION(e, msg);
            }
            return val;
        }

// ¸��Ƚ����ͤμ����ǽ�����2�ټ�֤ˤʤäƤ��뤬��
// ���ޡ��Ȥ˽���ˡ���פ��Ĥ��ʤ���
        bool Property::exist_(const string& key) const
        {
            if(dataMap_.find(key) == dataMap_.end() )
                return false;
            return true;
        }

// template��������
        template string Property::get(const string&) const;
        template string Property::get(const string&, const string&) const;
        template int    Property::get(const string&) const;
        template int    Property::get(const string&, const int&) const;
        template double Property::get(const string&) const;
        template double Property::get(const string&, const double&) const;
        template bool   Property::get(const string&) const;
        template bool   Property::get(const string&, const bool&) const;
        template char   Property::get(const string&) const;
        template char   Property::get(const string&, const char&) const;
// �ʤ������顼�Ȥʤ�ΤǼ�������������
//        template UbVectorD   Property::get(const string&) const;
//        template UbVectorD   Property::get(const string&, const char&) const;


    } // util
}

