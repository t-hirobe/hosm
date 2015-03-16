#ifndef EXCEPTION_EXCEPTION_H__
#define EXCEPTION_EXCEPTION_H__

/*! @brief
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 * Lastupdate: 2015-03-16 11:59:40
 */

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

/*! @brief すべての例外を補足
 *
 * このプロジェクトではこのクラスをすべての例外クラスの上位クラス
 * とし、すべての関数はこの派生クラスのみを送出します。
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 */
        class All
        {
        public:
            All() {};
            // エラーメッセージを追加
            void addInfo(const std::string& fileInfo, const std::string& msg);
            // 現在のエラーメッセージログを表示
            std::string info() const;

        private:
            std::vector<std::string> messages_;
        };


/*! @brief IO処理系の例外を補足
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 */
        class IO :
            public All
        {};


/*! @brief 解析系の例外を補足
 *
 * @author Tomoyuki Hirobe
 * @date 2010/05/03  5:17:37
 */
        class Parse :
            public All
        {};


    } // exception
} // flow

#endif /** UTIL_EXCEPTION_H__ */
