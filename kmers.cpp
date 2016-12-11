#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <sparsehash/sparse_hash_map>  //Google's sparse_hash_map is the most memory efficient
//using google::sparse_hash_map;
//#define _SECURE_SCL 0
int k,topKmers,HmapSize;
static void show_usage()
{
    std::cerr << "Usage: ./kmers --filename big.fastq --kmersize 30 --topcount 25\n"
              << "OR ./kmers -f big.fastq -k 30 -t 25\n"
			  << "Options:\n"
              //<< "--skipN or -sn : YES (Default is NO - skip Kmers that includes N)";
              << std::endl;
	//An option for Skipping K-mers including "N" is yet to be implemented: Reason N is I guess non-bases result
}

struct kMerStruct { //Vector placeholder
	std::string kmer;
	int count;
	bool operator() (const kMerStruct& str1 ,const kMerStruct& str2){ return(str1.count>str2.count);}
} kMerstr1; // object to make use of advantage of operators...

	//std::unordered_map< char *, int> KmersCounterMap(HmapSize);          //uses too much memory
	std::vector<kMerStruct> topKmersVector(topKmers+1);
	google::sparse_hash_map<std::string, int> KmersCounterMap(HmapSize);   // initizalization of the hash map with pre-calculated size

void initializer(std::string filename)
{
    //Hash maps and vectors are efficient if
    //they initialized with pre-determined sizes.
    //Therefore, here I retrieve the number of lines
    //predicate Hashmap size number according line and double it.
    std::cout<<"File "<<filename<<" is started processing with " <<k<<"-mers"<<std::endl;
    std::string line;
    std::ifstream fcounter;
    int countLines=1;
    fcounter.open(filename);
    while (fcounter.good())
    {
        getline(fcounter,line);
        countLines++;
    }
    fcounter.close();
    HmapSize=(((countLines/4)*(line.length()-k+1))*2);

}

void HmapToVector(std::string kmer, int count)
{
    int top=topKmers;
    bool KmerExist=false;
    kMerStruct KmerStr;
    KmerStr.kmer = kmer;
    KmerStr.count = count;
    for (int i = 0; i < topKmersVector.size(); i++) {
        if(kmer==topKmersVector.at(i).kmer){
            topKmersVector.at(i).count++;
            KmerExist=true;
            break;
        }
    }
    if(!KmerExist){
        topKmersVector.push_back(KmerStr);
    }
    // Vectors are time consuming, therefore
    // removing the last element after sorting is better
    // rather than sorting at the end of whole process.
    sort(topKmersVector.begin(), topKmersVector.end(), kMerstr1);
    if (topKmersVector.size() > top){
        topKmersVector.pop_back();
    }
}

void parseLinesToKmersViaHash(std::string  line)
{
	int kmerLength=k;
	std::string kMerString;
	for (int i = 0; i < line.length() - kmerLength +1; i++)  // +1 on condition is neccesary
	{
		kMerString=line.substr(i, kmerLength);
		if(KmersCounterMap.count(kMerString))
		{
			KmersCounterMap[kMerString]++;
		}
		else
		{
			KmersCounterMap.insert(std::pair<  std::string, int>(kMerString, 1));
		}
		HmapToVector(kMerString, KmersCounterMap[kMerString]);
	}
}

int main(int argc, char *argv[] )
{
	char* filename;
	std::string line;
	int counter = 3;
	show_usage();
	if (argc < 6) {
        show_usage();
        return 0;
    }
	std::ifstream file;
	//file.exceptions( std::ifstream::failbit | std::ifstream::badbit );

	for (int i = 1; i < argc; i++) {
       std::string arg= argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage();
            //return 0;
        } else if ((arg == "-f") || (arg == "--filename")) {
                filename = argv[++i];
                file.open(filename);
					if(!file.is_open()){std::cout<<"--filename is required" << std::endl; return 0;  }
               //try {} catch (const std::ifstream::failure& e) {}
        } else if ((arg == "-k") || (arg == "--kmersize")) {
                k = atoi(argv[++i]);
		} else if ((arg == "-t") || (arg == "--topcount")) {
                topKmers = atoi(argv[++i]);
		} 
	}
	initializer(filename);
	 while (file.good())
	{
		getline(file,line);

		if (counter % 4 == 0 && !line.empty() && !useVector){   //eof check at the end of file, do not leave the last line empty
			parseLinesToKmersViaHash(line);
		}
		counter++;
		if(counter%1000==0)
		{
			std::cout<<counter<<" lines are processed"<<std::endl;
		}
	}
	 file.close();
	 std::cout<<"Top 25 K-Mers"<<std::endl;
		for (int i = 0; i < topKmersVector.size() -1 ; i++)
		{
			kMerStruct topKmers = topKmersVector.at(i);
			std::cout << topKmers.kmer << ": " << topKmers.count << std::endl;
		}
	
	return 0;
}
