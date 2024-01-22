#include <string>
#include <cctype>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>


#include "alignmentIndex.h"
#include "mafParser.h"

using namespace std;

int main(int argc, char* argv[])
{
	std::string buffer;
	std::vector<std::string> query;
	/*
	std::string inFile("D:/Projects/JHU2/maf-query/data/subset.maf");
	std::string reference("Homo_sapiens");
	std::ifstream queryFile("D:/Projects/JHU2/maf-query/data/species.txt");
	*/

	std::string inFile(argv[1]);
	std::string reference(argv[2]);
	std::ifstream queryFile(argv[3]);

	while (std::getline(queryFile, buffer))
	{
		query.push_back(buffer);
	}

	std::sort(query.begin(), query.end());
	std::ifstream input(inFile.c_str());
	MafParser parser(input);
	std::vector<size_t> referenceIdx;
	std::vector<size_t> referencePos;
	std::vector<std::string> conservation;
	std::vector<MafParser::Record> alignmentBlock;

	AlignmentIndexReader index(query, "");
	while (parser.ReadBlock(alignmentBlock))
	{
		referenceIdx.clear();
		referencePos.clear();
		for (size_t i = 0; i < alignmentBlock.size(); ++i)
		{
			if (alignmentBlock[i].species == reference)
			{
				referenceIdx.push_back(i);
				referencePos.push_back(alignmentBlock[i].start);
			}
		}

		for (size_t j = 0; j < alignmentBlock[0].alignment.size(); j++)
		{
			for (size_t i = 0; i < referenceIdx.size(); ++i)
			{
				conservation.clear();
				auto& pos = referencePos[i];
				auto& ref = alignmentBlock[referenceIdx[i]];
				if (ref.alignment[j] == '-')
				{
					continue;
				}

				for (auto& rec : alignmentBlock)
				{
					if (rec.species != reference && toupper(ref.alignment[j]) == toupper(rec.alignment[j]))
					{
						conservation.push_back(rec.species);
					}
				}

				std::sort(conservation.begin(), conservation.end());
				conservation.erase(std::unique(conservation.begin(), conservation.end()), conservation.end());
				std::vector<size_t> indexConservation;
				index.Read(ref.chr, ref.PositivePos(pos), indexConservation);
				std::vector<std::string> strIndexConservation;
				for (auto idx : indexConservation)
				{
					strIndexConservation.push_back(query[idx]);
				}

				assert(strIndexConservation == conservation);
				++pos;
			}
		}
	}

	return 0;
}
