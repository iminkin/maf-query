#ifndef _MAF_PARSER_H_
#define _MAF_PARSER_H_

#include <string>
#include <istream>
#include <sstream>
#include <iostream>
#include <algorithm>

class MafParser
{
public:
	struct Record
	{
		std::string species;
		std::string chr;
		size_t start;
		size_t length;
		char strand;
		size_t chrSize;
		std::string alignment;

		size_t PositivePos(size_t pos)
		{
			if (strand == '+')
			{
				return pos;
			}

			return chrSize - pos - 1;
		}
	};

	MafParser(std::istream& input) : _input(input)
	{

	}

	bool ReadBlock(std::vector<Record>& retRecord)
	{		
		retRecord.clear();
		while (std::getline(_input, _buffer))
		{
			if (_buffer.size() != 0 && _buffer[0] == 'a')
			{
				while (std::getline(_input, _buffer))
				{
					if (_buffer.size() == 0)
					{
						return true;
					}
					else if (_buffer[0] == 's')
					{
						Record record;
						std::stringstream token(_buffer);
						token >> _buffer >> _buffer >> record.start >> record.length >> record.strand >> record.chrSize >> record.alignment;
						auto dot = std::find(_buffer.begin(), _buffer.end(), '.');
						record.species.assign(_buffer.begin(), dot);
						record.chr.assign(dot + 1, _buffer.end());
						retRecord.push_back(record);
					}
				}
			}			
		}

		return false;
	}

private:
	std::string _buffer;
	std::istream& _input;
};

#endif