/* Usage: ChangeCoordinate [-i <CGmap>] -I <Index>
 * Weilong GUO, start from 2016-03-09
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

struct parameter
{
	string infile;	// -i
	int chr_infile_col;	// -c
	int pos_infile_col;	// -p
	string indexfile;	// -I
	int chr_index_col;	// -C
	int pos_index_col;	// -P
	int target_index_col;	// -T
	string outfile;	// -o
};

parameter param;


void exit_with_help( void )
{
	printf(
		"Usage: ChangeCoordinate [-i <INPUT>] -I <INDEX>\n"
		"Author: Guo, Weilong; guoweilong@126.com; 2016-03-09\n"
		"Last updated : 2016-03-10\n"
		"Description: Changing the chr and pos according to the <Index> file.\n"
		"    This program is designed as succeeding program [ThreadFasta].\n"
		"    INPUT file: chr column can be in any order; pos column should be sorted.\n"
		"Options:\n"
		"  -i <STRING> : Input file, STDIN if omitted\n"
		"  -c <INT> : Column in Input file for chr name [Default: 1]\n"
		"  -p <INT> : Column in Input file for position [Default: 2]\n"
		"  -I <STRING> : Index file\n"
		"  -C <INT> : Column in Index file for original chr name [Default: 1]\n"
		"  -P <INT> : Column in Index file for original position [Default: 2]\n"
		"  -T <INT> : Column in Index file for target chr name [Default: 3]\n"
		//"  -o <STRING> : Output file\n"
		"  -h : help\n"
	);
	exit(1);
}

void parse_command_line(int argc, char **argv)
{
	int i;
	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 'i':
				param.infile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'c':
				param.chr_infile_col = atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'p':
				param.pos_infile_col = atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'I':
				param.indexfile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'C':
				param.chr_index_col = atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'P':
				param.pos_index_col = atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'T':
				param.target_index_col = atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			//case 'o':
			//	param.outfile = string(argv[i]);
			//	if(++i>argc) exit_with_help();
			//	break;
			case 'h':
				exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.indexfile.length() ){
		exit_with_help();
	}
}

class index_entry{
public:
	int pos;
	string target_chr;
};

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

void init(){
	// default values
	param.infile = "";
	param.chr_infile_col = 1;
	param.pos_infile_col = 2;
	param.indexfile	= "";
	param.chr_index_col	= 1;
	param.pos_index_col = 2;
	param.target_index_col = 3;
	param.outfile = "";
}


map<string, vector<index_entry> > INDEX;


int ReadIndexFile( string filename ) {
	int target_index_col = param.target_index_col-1;
	int chr_index_col = param.chr_index_col - 1;
	int pos_index_col = param.pos_index_col - 1;
	//cerr << "1\n";
	ifstream infile(filename.c_str());
	if(!infile) {
		cerr << "Error when reading input file : " << filename << ".\n";
		exit(-1);
	}
	string str;
	string pre_chr = "";
	while ( getline( infile, str ) ) {
		//cerr << "2\n";
		vector<string> tokens = string_tokenize(str);
		string chr = tokens[chr_index_col];
		int pos = atoi( tokens[pos_index_col].c_str() );
		string target_chr = tokens[target_index_col];
		index_entry NewEntry;
		NewEntry.pos = pos;
		NewEntry.target_chr = target_chr;
		//cerr << "3\n";
		//
		if ( pre_chr == chr ) {
			INDEX[chr].push_back(NewEntry);
		} else {
			map<string, vector<index_entry> >::iterator iter = INDEX.find(chr);
			if (iter == INDEX.end() ) { // New chr
				vector<index_entry> tmpVec;
				INDEX[chr] = tmpVec;
			}
			INDEX[chr].push_back(NewEntry);
			pre_chr = chr;
		}
	}
	infile.close();
	return 1;
}


int OutputIndexFile( void ) {
	//cerr << "5\n";
	map<string, vector<index_entry> >::iterator iter;
	for ( iter = INDEX.begin(); iter != INDEX.end(); iter++ ) {
		string chr = iter->first;
		vector<index_entry>::iterator viter;
		for ( viter = INDEX[chr].begin(); viter != INDEX[chr].end(); viter++ ) {
			cout << chr << "\t" << viter->pos << "\t" << viter->target_chr << endl;
		}
	}
	return 1;
}

string VectorJoin ( vector<string> & vec) {
	string str = "";
	if (vec.size() > 0) {
		vector<string>::iterator iter = vec.begin();
		str += (*iter);
		iter++;
		while ( iter != vec.end() ) {
			str += ("\t" + (*iter));
			iter++;
		}		
	}
	return str;
}

int ChangeCoordinate ( string filename )
{
	int chr_infile_col = param.chr_infile_col - 1;
	int pos_infile_col = param.pos_infile_col - 1;
	// Infile
	istream *p_infile = &cin;
	ifstream infile;
	if (filename != "") {
		infile.open(filename.c_str());
		if(!infile) {
			cerr << "cannot open input file" << filename.c_str() << endl;
			return -1;
		}
		p_infile = &infile;
	}
	/*ifstream infile(filename.c_str());
	if(!infile) {
		cerr << "Error when reading input file : " << filename << ".\n";
		return -1;
	}*/
	/*
	// Outfile
	ofstream of(param.outfile.c_str());
	if( !of ) {
		cerr << "Open output file error:" << param.outfile << "\n";
		exit(-1);
	}*/
	string str;
	string pre_chr = "";
	int pre_pos = 0;
	string pre_target = "";
	vector<index_entry>::iterator viter;
	// Reading file
	while ( getline( *p_infile, str ) ) {
		vector<string> tokens =  string_tokenize(str);
		string chr = tokens[chr_infile_col];
		int pos = atoi( tokens[pos_infile_col].c_str() );
		if (chr != pre_chr) {
			viter = INDEX[chr].begin();
			pre_pos = 0;
			pre_chr = chr;
		}
		if ( (viter != INDEX[chr].end()) && (pos >= viter->pos) ) {
			while ( (viter != INDEX[chr].end()) && (pos >= viter->pos) ) {
				pre_target = viter->target_chr;
				pre_pos = viter->pos;
				viter++;
			}
		}
		tokens[chr_infile_col] = pre_target;
		char buffer[100];
		sprintf(buffer, "%d", pos - pre_pos);
		tokens[pos_infile_col] = string(buffer);
		string msg = VectorJoin(tokens);
		cout << msg << endl;
	}
	//
	infile.close();
	cerr << "Finished reading FASTA file\n";
	return 1;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = "";
	param.chr_infile_col = 1;
	param.pos_infile_col = 2;
	param.indexfile	= "";
	param.chr_index_col	= 1;
	param.pos_index_col = 2;
	param.target_index_col = 3;
	param.outfile = "";
#else
	parse_command_line(argc,argv);
#endif

	ReadIndexFile( param.indexfile );
	//OutputIndexFile( );

	/*
	// test VectorJoin()
	vector<string> StrVec;
	StrVec.push_back("Guo");
	StrVec.push_back("Wei");
	StrVec.push_back("Long");
	cout << VectorJoin(StrVec) << endl;
	*/

	ChangeCoordinate( param.infile );


	return 1;
}

