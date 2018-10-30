/* Usage: ThreadFasta [-i <input.fa>]  -o <output.fa> -d <dict.bed>
 * Weilong GUO, start from 2016-03-03
 * *********************
 * Modification log:
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;
#include "math.h"

struct parameter
{
	string infile;	// -i
	string outfile;	// -o
	string dictfile;	// -d
	string chrname; // -c
	int MaxChrLen;	// -M
	int GapLen;	// -g
	bool Info;	// -I
};

class Seq{
public:
	string name;
	string content;
	Seq( void ) {
		name = "";
		content = "";
	}
	Seq(string NAME, string CONTENT) {
		name = NAME;
		content = CONTENT;
	}
};

parameter param;

void exit_with_help( void )
{
	printf(
		"Usage: ThreadFasta [-i <input.fa>]  -o <output.fa> -d <dict.bed>\n"
		"                   [-M 10000000 -g 50 -I]\n"
		"Author: Guo, Weilong; guoweilong@126.com; 2016-03-03\n"
		"Last updated: 2016-09-22\n"
		"Description: Threading short contigs into fewer large pseudo-chromosomes\n"
		"Options:\n"
		"  -i <STRING> : Input file, Fasta, STDIN if omitted\n"
		"  -o <STRING> : Output file, Fasta, constructed new file\n"
		"  -d <STRING> : Output file, TXT, dictionary for contig position\n"
		"  -c <STRING> : chromosone name group. [Default: None].\n" 
		"                Ex: AA, then chrAA_1, chrAA_2 will be used.\n"
		"  -g <INT> : Number of bases added between the contigs [Default: 50]\n"
		"  -M <INT> : maximum bp of constructed fasta sequence [Default: 10000000]\n"
		"  -I : Output interactive information when specified\n"
		"  -h : help\n"
	);
	exit(1);
}

vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
inline vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty) {
	// Skip delimiters at beginning.
	string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	vector<string> result;
	result.clear();
	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		//__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");
		//if (pos == lastPos) result.push_back("");
		result.push_back(str.substr(lastPos, pos - lastPos));
		if (pos == string::npos) break;
		if (pos == str.length() - 1) {
			if (!skip_empty) result.push_back("");
			break;
		}
		// Skip delimiters.  Note the "not_of"
		lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
	return result;
}

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);
}

void ToLowerString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))tolower);
}


map<string, string> chr2seq;

// todo: add gap with certain length

/** Output the structure in a FASTA-format file. By default, the length of each
 *  line is 50 bp.
 */
int OutputFasta( ofstream & of, Seq & seq, int step = 50 )
{
	if(!of) {
		cerr << "Error with the fasta output handle.\n";
		return -1;
	}
	of << ">" << seq.name << endl;
	int read_len = seq.content.length();
	int max_len = read_len - step;
	int i;
	//cout << max_len << endl;
	for (i=0; i<max_len; i+=step) {
		of << seq.content.substr(i, step) << endl;
	}
	if ( i<=(read_len-1) ) {
		of << seq.content.substr(i, read_len-i) << endl;
	}
	/*while( read.length() > len ) {
		of << read.substr(0,len) << endl;
		read = read.substr(len);
	}
	of << read << endl; */
	return seq.content.length();
}


int ThreadFasta ( string filename, vector<Seq> & SeqVec )
{
	// Initiation
	Seq contig_seq = Seq();
	Seq chr_seq = Seq();
	int chr_len = 0;
	int chr_id = 0;
	int chr_pos = 0;
	char buffer[100];
	string gap = "";
	//string gapchar = "N";
	//gap.append(gapchar.c_str(), param.GapLen);
	gap.append(param.GapLen, 'N');
	string str;
    cout << param.GapLen << "\t" << gap << endl;
	// Infile
	istream *p_infile = &cin;
	ifstream infile;
	if (filename != "") {
		infile.open(filename.c_str());
		if(!infile) {
			cerr << "cannot open input file : " << filename << endl;
			return -1;
		}
		p_infile = &infile;
	}
	// Outfile
	ofstream dict_of(param.dictfile.c_str());
	if( !dict_of ) {
		cerr << "Open output file error:" << param.dictfile << "\n";
		exit(-1);
	}
	ofstream fasta_of(param.outfile.c_str());
	if( !fasta_of ) {
		cerr << "Open output file error:" << param.outfile << "\n";
		exit(-1);
	}
	// Reading file
	if (param.Info) { cerr << "Reading FASTA file...\n"; }
	while ( getline( *p_infile, str ) ) {
		if( str[0] == '>' ) { // Reading the title line
			// For the First line of input file
			if ( chr_seq.name.empty() ) {
				// contig_seq
				vector<string> tokens = string_tokenize(str, "\t");
				contig_seq.name = tokens[0].substr(1, str.length());
				contig_seq.content = "";
				// chr_seq
				chr_id++; sprintf(buffer, "%d", chr_id);
				chr_seq.name = string("chr") + param.chrname + string("_") + buffer;
				continue;
			}
			int contig_len = contig_seq.content.length();
			if ( (chr_len + contig_len) >= param.MaxChrLen ) { // create new chr_seq
				// SeqVec.pushback(chr_seq);
				OutputFasta( fasta_of, chr_seq, 60);
				if (param.Info) { 
					cerr << "Finished " << chr_seq.name << ":\t";
					cerr << chr_seq.content.length() << endl; 
				}
				chr_seq = Seq();
				chr_id++; sprintf(buffer, "%d", chr_id);
				chr_seq.name = string("chr") + param.chrname + string("_")  + string(buffer);
			}
			// Attach contig to chr
			chr_pos = chr_seq.content.length() + 1;
			chr_seq.content += (contig_seq.content + gap);
			chr_len = chr_seq.content.length();
			// Store the contig position
			dict_of << chr_seq.name << "\t" << chr_pos << "\t" ;
			dict_of << contig_seq.name << endl;
			// Create new contig
			contig_seq = Seq();
			vector<string> tokens = string_tokenize(str, "\t");
			contig_seq.name = tokens[0].substr(1, str.length());
			contig_seq.content = "";
		} else { // Reading the seq line
			// Attach contig to chr
			contig_seq.content += str;
		}
	}
	// The last contig
	chr_pos = chr_seq.content.length() + 1;
	chr_seq.content += contig_seq.content;
	dict_of << chr_seq.name << "\t" << chr_pos << "\t" ;
	dict_of << contig_seq.name << endl;
	OutputFasta( fasta_of, chr_seq, 60);
	if (param.Info) { 
		cerr << "Finished " << chr_seq.name << ":\t";
		cerr << chr_seq.content.length() << " bp" << endl; 
	}
	// End
	dict_of.close();
	fasta_of.close();
	infile.close();
	cout << "Finished reading FASTA file\n";
	return 1;
}

void parse_command_line(int argc, char **argv)
{
	int i;
	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'i':
				param.infile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'd':
				param.dictfile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'c':
				param.chrname = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'M':
				param.MaxChrLen= atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'g':
				param.GapLen= atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'I':
				param.Info = true;
				break;
			case 'h':
				exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.outfile.length() || !param.dictfile.length() ){
		exit_with_help();
	}
}



void init(){
	// default values
	param.infile = "";
	param.outfile = "";
	param.dictfile = "";
	param.chrname = "";
	param.MaxChrLen	= 1000000;
	param.GapLen	= 50;
	param.Info = false;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.outfile = string("output.txt");
	param.MaxChrLen	= 1000000;
	param.GapLen	= 50;
#else
	parse_command_line(argc,argv);
#endif
	//GetGenome( param.infile.c_str() );

	vector<Seq> SeqVec;	// Vector to store the sequences
	ThreadFasta( param.infile, SeqVec );


	return 1;
}

