#ifndef UTIL_PROPERTY_H__
#define UTIL_PROPERTY_H__

/*! @brief 設定ファイル簡易読み込みクラス
 * @author Tomoyuki Hirobe
 * @date 2010/05/02 19:56:08
 * Lastupdate: 2015-03-16 11:55:42
 */

#include <iostream>
#include <map>
#include <string>

namespace flow {
namespace util {


/*! @brief 設定ファイル簡易読み込みクラス
 *
 * "="によって分割された左辺値をキー、右辺値を値として
 * ファイルから読みこむ。<br/>
 * "#"で始まる行はコメント行として読み飛ばす
 * @author Tomoyuki Hirobe
 * @date 2010/05/02 19:56:08
 */
class Property
{
public:
    /// デフォルトコンストラクタ
    explicit Property();

    /*! @brief 初期化時に指定ファイル名のファイルを設定ファイルとして読み込む。
     * @param fileName 読み込み設定ファイル名
     * @throw IOException ファイル読み込み失敗時
     */
    explicit Property(const std::string& fileName);

    /// デストラクタ
    virtual ~Property();

    /*! @brief 指定ファイル名のファイルを設定ファイルとして読み込む。
     * @param fileName 読み込み設定ファイル名
     * @throw IOException ファイル読み込み失敗時
     */
    void read(const std::string& fileName);


    /*! @brief 入力キー文字列に対応する値をType型で取得する。
     *
     * @param key キー文字列
     * @return キーに対応する値
     * @throw IOException キーが存在しない場合
     * @throw ParseException 対応値がType型に変換できない場合
     */
    template<typename Type>
    Type get(const std::string& key) const;


    /*! @brief 入力キー文字列に対応する値をType型で取得する。
     *
     * @param key キー文字列
     * @param defaultVal キーが存在しない場合にこの値を返す
     * @return キーに対応する値、キーが存在しない場合はdefaultVal
     * @throw ParseException 対応値がType型に変換できない場合
     */
    template<typename Type>
    Type get(const std::string& key, const Type& defaultVal) const;

    /// デバッグ用、読み込んだ全情報をエラー出力に表示
    void output() const;

private:
    std::string fileName_;
    std::map<std::string, std::string> dataMap_;
    bool exist_(const std::string& key) const;
    void setData_(std::ifstream& in);
    template<typename Type>
    Type cast_(const std::string& valStr) const;
};


} // util
} // flow

#endif /** UTIL_PROPERTY_H__ */
