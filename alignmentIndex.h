#ifndef _ALIGNMENT_INDEX_H_
#define _ALIGNMENT_INDEX_H_


#include <map>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <algorithm>

class AlignmentIndexWriter
{
public:
	AlignmentIndexWriter(const std::vector<std::string>& genomes, const std::string& reference, const std::string & prefix):
		_genomes(genomes),
		_reference(reference),
		_prefix(prefix),
		_blockSize(0)
	{
		auto bits = sizeof(BlockT) * 8;
		_blockSize = genomes.size() / bits;
		if (genomes.size() % bits)
		{
			++_blockSize;
		}

		_buffer.reset(new BlockT[_blockSize]);
		
	}

	void Write(const std::string& chr, size_t pos, const std::vector<std::string>& conservation)
	{
		auto it = _file.find(chr);
		if (it == _file.end())
		{			
			_file[chr].reset(new std::ofstream(_prefix + chr, std::ios::binary | std::ios::out));
			it = _file.find(chr);
		}

		size_t bit;
		size_t element;
		auto ptr = _buffer.get();
		std::fill(ptr, ptr + _blockSize, 0);
		for (auto& g : conservation)
		{
			size_t idx = std::lower_bound(_genomes.begin(), _genomes.end(), g) - _genomes.begin();
			GetElementBit(idx, element, bit);
			_buffer[element] |= (BlockT(1) << bit);	
		}

		size_t block = pos * _blockSize * sizeof(BlockT);
		it->second->seekp(block);
		it->second->write(reinterpret_cast<const char*>(ptr), sizeof(_buffer[0]) * _blockSize);
	}

private:
	
	void GetElementBit(size_t idx, size_t& element, size_t & bitPos)
	{
		auto bits = sizeof(BlockT) * 8;
		element = idx / bits;
		bitPos = idx % bits;
	}

	std::vector<std::string> _genomes;
	std::string _reference;
	std::string _prefix;
	typedef std::unique_ptr<std::ofstream> File;
	std::map<std::string, File> _file;	
	typedef uint32_t BlockT;
	std::unique_ptr<BlockT[]> _buffer;
	size_t _blockSize;
};


class AlignmentIndexReader
{
public:
	AlignmentIndexReader(const std::vector<std::string>& genomes, const std::string& prefix) :
		_prefix(prefix),
		_blockSize(0)
	{
		auto bits = sizeof(BlockT) * 8;
		_blockSize = genomes.size() / bits;
		if (genomes.size() / bits)
		{
			++_blockSize;
		}

		_buffer.reset(new BlockT[_blockSize]);

	}

	void Read(const std::string& chr, size_t pos, std::vector<size_t>& conservation)
	{
		auto it = _file.find(chr);
		if (it == _file.end())
		{
			_file[chr].reset(new std::ifstream(_prefix + chr, std::ios::binary | std::ios::in));
			it = _file.find(chr);
		}

		conservation.clear();
		auto ptr = _buffer.get();
		size_t block = pos * _blockSize * sizeof(BlockT);
		std::fill(ptr, ptr + _blockSize, 0);
		it->second->seekg(block);
		it->second->read(reinterpret_cast<char*>(ptr), sizeof(_buffer[0]) * _blockSize);
		for (size_t i = 0; i < _blockSize; i++)
		{
			size_t bits = sizeof(BlockT) * 8;
			for (size_t j = 0; j < bits; j++)
			{				
				if (_buffer[i] & (BlockT(1) << j))
				{
					conservation.push_back(bits * i + j);
				}
			}
		}
	}

private:
	std::string _prefix;
	typedef std::unique_ptr<std::ifstream> File;
	std::map<std::string, File> _file;
	typedef uint32_t BlockT;
	std::unique_ptr<BlockT[]> _buffer;
	size_t _blockSize;
};

#endif
