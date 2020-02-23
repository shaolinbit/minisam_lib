#ifndef __CSV_PARSER_H__
#define __CSV_PARSER_H__

#include <ios>
#include <vector>
#include <string>


// 当前直接实现成将数据解析成double

class CSVParser
{

    typedef std::vector<double> ROW_TYPE;
    typedef std::vector<ROW_TYPE*> DATA_TYPE;

public:
    CSVParser() {};
    virtual ~CSVParser();

    // return rows，强制认为第一个非空行就是标题行
    size_t load(const std::string &csvFile,
                const std::string &sep,
                const std::vector<std::string> &header);

    size_t size() const
    {
        return data.size();
    }

    const ROW_TYPE& getLine(int i) const
    {
        return *data.at(i);
    }
    DATA_TYPE::const_iterator begin() const
    {
        return data.begin();
    }
    DATA_TYPE::const_iterator end() const
    {
        return data.end();
    }

    std::string& trim(std::string& s, const std::string& char_list)
    {
        if (!s.empty())
        {
            s.erase(0, s.find_first_not_of(char_list));
            s.erase(s.find_last_not_of(char_list) + 1);
        }
        return s;
    }
    std::vector<std::string> &split(const std::string &s,
                                    std::vector<std::string> &ret,
                                    const std::string &sep);

private:
    DATA_TYPE data;
};


#endif
