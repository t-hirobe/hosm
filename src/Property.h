#ifndef UTIL_PROPERTY_H__
#define UTIL_PROPERTY_H__

#include <iostream>
#include <map>
#include <string>

namespace flow {
namespace util {


/*! @brief ����ե�����ʰ��ɤ߹��ߥ��饹
 *
 * "="�ˤ�ä�ʬ�䤵�줿�����ͤ򥭡��������ͤ��ͤȤ���
 * �ե����뤫���ɤߤ��ࡣ<br/>
 * "#"�ǻϤޤ�Ԥϥ����ȹԤȤ����ɤ����Ф�
 * @author Tomoyuki Hirobe
 * @date 2010/05/02 19:56:08
 * $Lastupdate: 2010/05/16 23:43:07 $
 */
class Property
{
public:
    /// �ǥե���ȥ��󥹥ȥ饯��
    explicit Property();

    /*! @brief ��������˻���ե�����̾�Υե����������ե�����Ȥ����ɤ߹��ࡣ
     * @param fileName �ɤ߹�������ե�����̾
     * @throw IOException �ե������ɤ߹��߼��Ի�
     */
    explicit Property(const std::string& fileName);

    /// �ǥ��ȥ饯��
    virtual ~Property();

    /*! @brief ����ե�����̾�Υե����������ե�����Ȥ����ɤ߹��ࡣ
     * @param fileName �ɤ߹�������ե�����̾
     * @throw IOException �ե������ɤ߹��߼��Ի�
     */
    void read(const std::string& fileName);


    /*! @brief ���ϥ���ʸ������б������ͤ�Type���Ǽ������롣
     *
     * @param key ����ʸ����
     * @return �������б�������
     * @throw IOException ������¸�ߤ��ʤ����
     * @throw ParseException �б��ͤ�Type�����Ѵ��Ǥ��ʤ����
     */
    template<typename Type>
    Type get(const std::string& key) const;


    /*! @brief ���ϥ���ʸ������б������ͤ�Type���Ǽ������롣
     *
     * @param key ����ʸ����
     * @param defaultVal ������¸�ߤ��ʤ����ˤ����ͤ��֤�
     * @return �������б������͡�������¸�ߤ��ʤ�����defaultVal
     * @throw ParseException �б��ͤ�Type�����Ѵ��Ǥ��ʤ����
     */
    template<typename Type>
    Type get(const std::string& key, const Type& defaultVal) const;

    /// �ǥХå��ѡ��ɤ߹����������򥨥顼���Ϥ�ɽ��
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

