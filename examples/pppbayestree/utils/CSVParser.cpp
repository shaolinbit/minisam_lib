#include "./CSVParser.h"

#include <fstream>
#include <sstream>

using namespace std;

CSVParser::~CSVParser()
{
    for (DATA_TYPE::iterator it=data.begin(); it!=data.end(); ++it)
    {
        delete *it;
    }
    data.clear();
}

size_t CSVParser::load(const string &csvFile,
                       const string &sep,
                       const vector<string> &header)
{
    size_t rows = 0;
    ifstream fin(csvFile);
    string line;
    string char_list(" \n\r\t");
    vector<string> items;

    while (getline(fin, line))
    {
        trim(line, char_list);
        // 空行
        if (line.empty())
        {
            continue;
        }
        items.clear();
        split(line, items, sep);
        // 不规则的行
        if (header.size() > 0 && items.size() != header.size())
        {
            continue;
        }
        // 判断是不是 header
        if (rows++ == 0 && header.size() > 0 && items[0] == header[0])
        {
            continue;
        }
        // 转换数据
        ROW_TYPE *row_data = new ROW_TYPE;
        for (vector<string>::const_iterator it=items.begin(); it!=items.end(); ++it)
        {
            istringstream iss(*it);
            double tmp = 0.0;
            iss >> tmp;
            row_data->push_back(tmp);
        }
        data.push_back(row_data);
    }
    return rows;
}

vector<string>& CSVParser::split(const string& s, vector<string>& ret, const string& sep)
{
    size_t start = 0;
    size_t end = s.find_first_of(sep);
    while (end != string::npos)
    {
        ret.push_back(s.substr(start, end - start));
        start = end + 1;
        end = s.find_first_of(sep, start);
    }
    if (end > start)
    {
        ret.push_back(s.substr(start, end - start));
    }
    return ret;
}
