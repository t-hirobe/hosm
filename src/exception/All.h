#ifndef EXCEPTION_EXCEPTION_H__
#define EXCEPTION_EXCEPTION_H__

#include <string>
#include <vector>
#include <boost/format.hpp>
#include <boost/current_function.hpp>

#define MACRO_FILE_INFO                                     \
  (boost::format("%s(%d): %s")                              \
   % __FILE__ % __LINE__ % std::string(BOOST_CURRENT_FUNCTION)).str()    \

#define MACRO_THROW_EXCEPTION(e, msg)        \
  e.addInfo( MACRO_FILE_INFO, msg);          \
  throw


namespace flow {
    namespace exception {

/*! @brief ���٤Ƥ��㳰����­
 *
 * ���Υץ������ȤǤϤ��Υ��饹�򤹤٤Ƥ��㳰���饹�ξ�̥��饹
 * �Ȥ������٤Ƥδؿ��Ϥ����������饹�Τߤ����Ф��ޤ���
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 * $Lastupdate: 2010/05/06  2:17:26 $
 */
        class All
        {
        public:
            All() {};
            // ���顼��å��������ɲ�
            void addInfo(const std::string& fileInfo, const std::string& msg);
            // ���ߤΥ��顼��å���������ɽ��
            std::string info() const;

        private:
            std::vector<std::string> messages_;
        };


/*! @brief IO�����Ϥ��㳰����­
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 * $Lastupdate: 2010/05/06  2:17:26 $
 */
        class IO :
            public All
        {};


/*! @brief ���ϷϤ��㳰����­
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 * $Lastupdate: 2010/05/06  2:17:26 $
 */
        class Parse :
            public All
        {};


    } // exception
} // flow

#endif /** EXCEPTION_EXCEPTION_H__ */

