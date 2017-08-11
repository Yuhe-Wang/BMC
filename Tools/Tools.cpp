#include "Tools.h"
#include "lz4.h"
#include "quicklz.h"

bool existFile(const char* filename)
{
	FILE* fp = fopen(filename, "rb");
	if (fp)
	{
		fclose(fp);
		return true;
	}
	return false;
}
long fileSize(FILE *fp)
{
	long save_pos;
	long size_of_file;

	/* Save the current position. */
	save_pos = ftell(fp);

	/* Jump to the end of the file. */
	fseek(fp, 0L, SEEK_END);

	/* Get the end position. */
	size_of_file = ftell(fp);

	/* Jump back to the original position. */
	fseek(fp, save_pos, SEEK_SET);

	return size_of_file;
}
void replaceAll(string& str, const string& from, const string& to)
{
	if (from.empty())
		return;
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos)
	{
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
	}
}
int get_thread_num()
{
#ifdef USE_OPENMP
#ifdef WIN32 
	SYSTEM_INFO info;
	GetSystemInfo(&info);
	return info.dwNumberOfProcessors;
#else //linux/unix
	return get_nprocs();
#endif
#else
	return 1;
#endif
}
double nodeSpeedTest(int NThread) //return seconds to finish this function
{
	//processWeight = nodeSpeedTest(NThread);
	RunTimeCounter nodeRC;
	//do some parallel heavy job
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NThread)
#endif
	{
		double x = 0;
		for (volatile double i = 0; i < 1e7; i += 1.0)
		{
			x += sin(i);
			x -= cos(i);
		}
		x++;
	}
	return nodeRC.stop();
}
void exitApp(const char *inf)
{
	printf("fatal error: %s", inf);
	printf("\nPress enter key to exit...");
	getchar();
	exit(-1);
}
void matchFile(const string& dir, string& match)
{
	fs::recursive_directory_iterator end_itr;
	string filter = dir;
	filter += FileSeparator;
	filter += match;
	for (fs::recursive_directory_iterator ifile(dir); ifile != end_itr; ++ifile)
	{
		// Skip if not a file
		if (!fs::is_regular_file(ifile->status())) continue;

		if (ifile->path().string().find(filter) != std::string::npos)
		{
			match = ifile->path().string();
			return; //only find the first match
		}
	}
	match = ""; //means no match
}

string changeExtension(const char* fname, const char* newExt)
{
	string str(fname);
	size_t idot = str.find_last_of(".");
	if (idot == string::npos) return string();
	str = str.substr(0, idot+1);
	if (newExt[0] == '.') str += (newExt + 1);
	else  str += newExt;
	return str;
}

string findFirstFileInDir(const string& dir, const string& match) 
{
	// Match must use regex, i.e. you should use .* instead of *, and \\. instead of .
	// and .? instead of ?
	// For instance, you should use ".*\\.txt" instead of "*.txt"
	regex pattern(match);
	fs::directory_iterator end_itr;
	for (fs::directory_iterator jfile(dir); jfile != end_itr; ++jfile)
	{
		if (!fs::is_regular_file(jfile->status())) continue; // Skip if not a file
		// only get the file name
		if (std::regex_match(jfile->path().filename().string(), pattern))
		{
			return jfile->path().string();
		}
	}
	return string();
}

void pauseSeconds(const unsigned int t_seconds)
{
#ifdef WIN32
	Sleep(t_seconds * 1000);
#else //linux
	sleep(t_seconds);
#endif // WIN32
}

string printWithComma(const double fnum)
{
	char ch[128];
	int NExp = int(log(fnum) / log(10));
	if (NExp < 0) return string();
	char* buff = ch;
	double f = fnum;//change it to non-constant
	for (int i = 0; i <= NExp; ++i)
	{
		int dgt = int(f / exp((NExp - i)*log(10)));
		sprintf(buff, "%d", dgt);
		++buff;
		f -= dgt*exp((NExp - i)*log(10));
		if ((NExp - i) % 3 == 0 && i != NExp)
		{
			sprintf(buff, ",");
			++buff;
		}
	}
	return string(ch);
	
}
/******************************* Begin: LogProvider *************************************/
LogProvider Log;
std::mutex logMutex;
LogProvider::LogProvider()
{
	fp = NULL;
	b_shortMode = false;
	callBackFunc = NULL;
}
LogProvider::~LogProvider()
{
	if (fp) fclose(fp);
	//if (pBuffer) delete[] pBuffer;
}
void LogProvider::setLogName(int i, const string mode, string cdir) // the log name would be Process-i.log
{
	char name[60];
	sprintf(name, "Process-%d.log", i);
	cdir += name;
	if (fp != NULL) fclose(fp);
	fp = fopen(cdir.c_str(), mode.c_str());
	logName = cdir;

	gethostname(name, 60);
	printf("The host name = %s\n\n", name);
	fprintf(fp, "The host name = %s\n\n", name);
}
void LogProvider::setLogName(const char* name)
{
	logName = name;
	if (fp != NULL) fclose(fp);
	fp = fopen(name, "w");
	char hname[60];
	gethostname(hname, 60);
	printf("The host name = %s\n\n", hname);
	fprintf(fp, "The host name = %s\n\n", hname);
}
void LogProvider::closeFile()
{
	if (fp) fclose(fp);
	fp = NULL;
}
void LogProvider::flush()
{
	fflush(stdout);
	if (fp) fflush(fp);
}
void LogProvider::operator() (const char *fmt, ...)
{
	logMutex.lock();
	char logBuffer[4096];
	va_list params;
	va_start(params, fmt);
	int end = vsprintf(logBuffer, fmt, params);
	if (end > 0)
	{
		end--;
		if (logBuffer[end] == '\n')
		{
			logBuffer[end] = '\0';
		}
	}
	else
	{
		sprintf(logBuffer, "InformationLog: Invalid Information log parameters");
	}

	if (callBackFunc) callBackFunc(logBuffer); // it will replace the printing to screen
	else
	{
		if (b_shortMode)
		{
			printf("\r                                                                          ");
			printf("\r%s", logBuffer);
		}
		else printf("%s\n", logBuffer); //output and point to next line
		fflush(stdout);
	}
	if (fp)
	{
		fprintf(fp, "%s\n", logBuffer);
		fflush(fp);
	}
	va_end(params);
	logMutex.unlock();
}
char* LogProvider::timeNow()
{
	time_t nowtime = time(NULL);
	return ctime(&nowtime);
}

void LogProvider::warningColor() //change the text to alerting color
{
#ifdef WIN32
	HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	if (handle == 0) return;
	SetConsoleTextAttribute(handle, FOREGROUND_RED | FOREGROUND_BLUE);
#endif
}
void LogProvider::normalColor() //change the text to normal color
{
#ifdef WIN32
	HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	if (handle == 0) return;
	SetConsoleTextAttribute(handle, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE); //white color
#endif
}

void LogProvider::shortMode(bool mode)// only output in one line on screen 
{
	if (b_shortMode != mode && mode == false) printf("\n"); //need to jump to next line when switching from short to long mode
	b_shortMode = mode;
}

void LogProvider::setCallBack(LogCallBack pfunc)
{
	callBackFunc = pfunc;
}
/*******************************  End:  LogProvider *************************************/

/******************************* Begin: PRNG ********************************************/
#include "PRNG_CPP.h" //include the implementation here
/*******************************  End:  PRNG ********************************************/


/******************************* Begin: ConfigFile **************************************/
ConfigFile::ConfigFile(const string &input) { iLast = -1; if (!parse(input)) exitApp("fail to parse the config file!"); }

ConfigFile::ConfigFile(const char *filename) { iLast = -1; if (!parse(filename)) exitApp("fail to parse the config file!"); }

ConfigFile::~ConfigFile()
{
	for (unsigned int it = 0; it != blocks.size(); ++it)
		delete blocks[it].second;
}
ConfigFile* ConfigFile::getBlock(const string &key, bool followLast)
{
	for (unsigned int i = followLast ? iLast + 1 : 0; i < blocks.size(); ++i)
	{
		if (0 == key.compare(blocks[i].first)) //found the key
		{
			iLast = i;
			return blocks[i].second;
		}
	}
	return NULL;
}
vector<ConfigFile*> ConfigFile::getBlockArray(const string &key)
{
	vector<ConfigFile*> res;
	for (unsigned int i = 0; i < blocks.size(); ++i)
	{
		if (0 == key.compare(blocks[i].first)) //found the key
		{
			res.push_back(blocks[i].second);
		}
	}
	return res;
}
bool ConfigFile::getValue(const string &key, string &value, bool followLast)
{
	for (unsigned int i = followLast ? iLast + 1 : 0; i < mps.size(); ++i)
	{
		if (0 == key.compare(mps[i].first)) //found the key
		{
			value = mps[i].second;
			iLast = i;
			return true;
		}
	}
	return false;
}
bool ConfigFile::getValue(const string &key, int &value, bool followLast)
{
	for (unsigned int i = followLast ? iLast + 1 : 0; i < mps.size(); ++i)
	{
		if (0 == key.compare(mps[i].first)) //found the key
		{
			//parse the value part
			stringstream in(mps[i].second);
			in >> value;
			if (!in.fail()) { iLast = i; return true; }
		}
	}
	return false;
}
bool ConfigFile::getValue(const string &key, double &value, bool followLast)
{
	for (unsigned int i = followLast ? iLast + 1 : 0; i < mps.size(); ++i)
	{
		if (0 == key.compare(mps[i].first)) //found the key
		{
			//parse the value part
			stringstream in(mps[i].second);
			in >> value;
			if (!in.fail())
			{
				iLast = i;
				return true;
			}
		}
	}
	return false;
}
bool ConfigFile::getValue(const string &key, vector<string> &values, bool followLast)
{
	values.resize(0);
	for (unsigned int i = followLast ? iLast + 1 : 0; i < mps.size(); ++i)
	{
		if (0 == key.compare(mps[i].first)) //found the key
		{
			//parse the value part
			string item;
			stringstream in(mps[i].second);

			while (in.good())
			{
				in >> item;
				if (!in.fail() && !item.empty()) values.push_back(item);
			}
			if (values.size() > 0) { iLast = i; return true; }
		}
	}
	return false;
}
bool ConfigFile::getValue(const string &key, vector<int> &values, bool followLast)
{
	values.resize(0);
	for (unsigned int i = followLast ? iLast + 1 : 0; i < mps.size(); ++i)
	{
		if (0 == key.compare(mps[i].first)) //found the key
		{
			//parse the value part
			int item;
			stringstream in(mps[i].second);

			while (true)
			{
				in >> item;
				if (!in.fail()) values.push_back(item);
				else break;
			}
			if (values.size() > 0) { iLast = i; return true; }
		}
	}
	return false;
}
bool ConfigFile::getValue(const string &key, vector<double> &values, bool followLast)
{
	values.resize(0);
	for (unsigned int i = followLast ? iLast + 1 : 0; i < mps.size(); ++i)
	{
		if (0 == key.compare(mps[i].first)) //found the key
		{
			//parse the value part
			double item;
			stringstream in(mps[i].second);

			while (true)
			{
				in >> item;
				if (!in.fail()) values.push_back(item);
				else break;
			}
			if (values.size() > 0) { iLast = i; return true; }
		}
	}
	return false;
}

bool ConfigFile::parse(const string &input) //parse the content
{
	stringstream in(input);
	char ch;
	string key, value;

	while (true)
	{
		//get first printable character 
		while (true)
		{
			in.get(ch);
			if (33 <= ch&&ch <= 126) break;
			if (in.eof() || in.fail()) return true;
		}
		if ('#' == ch || '!' == ch || '%' == ch) while (in.get(ch), ch != '\n'); //skip this comment line
		else //effective data, deal with n lines
		{
			key = ch;
			while (in.get(ch), ch != '=')
			{
				if (in.eof() || !in.good()) return false;
				key += ch;
			}
			//check  if the value is empty in this line
			value = "";
			while (in.get(ch), ch != '\n') //process the rest of this line
			{
				if ('#' == ch || '!' == ch || '%' == ch)
				{
					while (in.get(ch), ch != '\n'); //skip this line
					break;
				}
				else value += ch;
			};

			trim(key);
			trim(value);

			if (0 == value.length()) //it means there should follow a brace bracket defining a child block
			{
				while (in.get(ch), ch != '{');
				string subSection;
				int bracePair = 1;
				while (true)
				{
					in.get(ch);
					if ('{' == ch) ++bracePair;
					if ('}' == ch) --bracePair;
					if (bracePair != 0) subSection += ch;
					else break;
				}
				while ((in.get(ch), ch != '\n') && !in.eof()); //skip the rest of this line

				ConfigFile * cf = new ConfigFile(subSection);
				blocks.push_back(make_pair(key, cf));
			}
			else
			{
				mps.push_back(make_pair(key, value));
			}
		}
	}
}
bool ConfigFile::parse(const char *filename)  //open the file and parse the whole content
{
	FILE *fp = fopen(filename, "r");
	if (NULL == fp) return false;
	char ch;
	string content;
	do
	{
		ch = getc(fp);
		if ('#' == ch || '!' == ch || '%' == ch)
		{
			while (getc(fp) != '\n');//skip this comment line;
			content += '\n';
			continue;
		}
		content += ch;
	} while (ch != EOF);
	if (fp) fclose(fp);
	return parse(content);
}

void ConfigFile::macroReplace(const MPS &ps)
{
	for (unsigned int j = 0; j < ps.size(); ++j) //replace all keys
	{
		for (unsigned int i = 0; i < mps.size(); ++i) replaceAll(mps[i].second, ps[j].first, ps[j].second);
	}

	for (unsigned int i = 0; i < blocks.size(); ++i)
	{
		blocks[i].second->macroReplace(ps);
	}

}
void ConfigFile::macroReplace(const string &from, const string& to)
{
	for (unsigned int i = 0; i < mps.size(); ++i) replaceAll(mps[i].second, from, to);
	for (unsigned int i = 0; i < blocks.size(); ++i)
	{
		blocks[i].second->macroReplace(from, to);
	}
}

void ConfigFile::trim(string& s)
{
	if (0 == s.size()) return;
	int size = (int)s.length();
	string ss;
	int left = size, right = -1;
	for (int i = 0; i < size; ++i)
	{
		if (33 <= s[i] && s[i] <= 126)
		{
			left = i;
			break;
		}
	}
	for (int i = size - 1; i >= 0; --i)
	{
		if (33 <= s[i] && s[i] <= 126)
		{
			right = i;
			break;
		}
	}
	if (left <= right)
	{
		for (int i = left; i <= right; ++i) ss += s[i];
		s = ss;
	}
	else s = "";
}
/*******************************  End:  ConfigFile **************************************/

/******************************* Begin: EmailNotify *************************************/
void EmailNotify::init(ConfigFile* cf)
{
	string str;
	cf->getValue("email notification", str);
	if (str.compare("yes") == 0) notify = true;
	else notify = false;
	if (cf->getValue("email encrypted file", str) && getEncryptedAccount(str)) return;
	cf->getValue("email address", address);
	cf->getValue("email password", pwd);
}
bool EmailNotify::getEncryptedAccount(string& filename)
{
	char str[200];
	FILE* fp = fopen(filename.c_str(), "rb");
	if (NULL == fp) return false;
	fread(str, sizeof(char)*200, 1, fp);
	fclose(fp);

	//decryption
	const char* key = "simpleXORcipher";
	int keylen = sizeof(key) / sizeof(char);
	for (int i = 0; i < 200; ++i)
	{
		for (int j = keylen-1; j >=0; --j)
		{
			str[i] = str[i] ^ key[j];
		}
	}

	//break the string to account and password
	int ind = 0;
	while (str[ind] != 0) ++ind;
	address = str;
	pwd = str + ind + 1;

	return true;
}
bool EmailNotify::doNotify(){ return notify; }
void EmailNotify::send(bool normal, LogProvider* pLog)
{
	char charBuff[256];
	gethostname(charBuff, 256);
	string hostName = charBuff;

	//create email.txt first
	string content = "From: <";
	content += address;
	content += ">\nTo: <";
	content += address;
	if (normal) content += ">\nSubject: simulation finished on computer ";
	else content += ">\nSubject: simulation exited abnormally on computer ";
	content += hostName;
	if (normal) content += "!\n\nDear user, \n\nThe simulation task was finished on the computer named ";
	else content += "!\n\nDear user, \n\nThe simulation task exited abnormally on the computer named ";
	content += hostName;
	content += ". Please login to check the result.\n\nBest wishes!";
	if (pLog != NULL)
	{
		FILE* fp = fopen(pLog->getLogName().c_str(), "r");
		if (fp != NULL)
		{
			content += "\n\n---------------------------------- Last Log Content -------------------------------------\n\n";
			char ch = fgetc(fp);
			while (ch != EOF)
			{
				content += ch;
				ch = fgetc(fp);
			}
			fclose(fp);
		}
	}
	FILE* fp = fopen("mail.txt", "w");
	if (NULL == fp) exitApp("Cannot write the email to file!");
	fprintf(fp, "%s", content.c_str());
	fclose(fp);

	//send the email vi curl
	string cmd;
#ifdef WIN32
	cmd = "curl --url \"smtps://smtp.gmail.com:465\" --ssl-reqd --mail-from \"";
#else //linux
	cmd = "curl --url \"smtps://smtp.gmail.com:465\" --ssl-reqd --mail-from \"";
#endif // WIN32
	cmd += address;
	cmd += "\" --mail-rcpt \"";
	cmd += address;
	cmd += "\" --upload-file mail.txt --user \"";
	cmd += address;
	cmd += ":";
	cmd += pwd;
	cmd += "\" --insecure";
	system(cmd.c_str());
}
/*******************************  End:  EmailNotify *************************************/

/******************************* Begin: BinaryFile **************************************/

BinaryFile::~BinaryFile()
{
// 	for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
// 	{
// 		if (it->second.len > 0) delete[] it->second.p;
// 	}
}
bool BinaryFile::read(const char* fname, bool onlyTest)
{
	FILE* fp = fopen(fname, "rb");
	if (NULL == fp) return false;
	delAll(); //remove all previous data

	//read the header 
	if (!header.read(fp)) return false;
	if (onlyTest)
	{
		fclose(fp);
		return true;
	}

	int NMap = 0;
	char cname[64];
	Chunk ck;
	
	if (header.b_compress)
	{
		int32_t original_size;
		int32_t compressed_size;
		fread(&original_size, sizeof(int32_t), 1, fp);
		fread(&compressed_size, sizeof(int32_t), 1, fp);
		char* source = new char[compressed_size];
		fread(source, compressed_size, 1, fp);	
		char* destination = new char[original_size];
		//decompressing
		if (0 == header.compressMethod)
		{
			LZ4_decompress_fast(source, destination, original_size);
		}
		else if (1 ==header.compressMethod)
		{
			qlz_state_decompress *qlz_state = new qlz_state_decompress;
			qlz_decompress(source, destination, qlz_state);
			delete qlz_state;
		}
		//now begin reading data from char* destination
		size_t cur = 0;
		memcpy(&NMap, destination + cur, sizeof(int32_t));
		cur += sizeof(int32_t);
		for (int i = 0; i < NMap; ++i)
		{
			memcpy(cname, destination + cur, 64 * sizeof(char));
			cur += 64 * sizeof(char);
			string tname((char*)cname);
			int32_t len;
			memcpy(&len, destination + cur, sizeof(int32_t));
			cur += sizeof(int32_t);
			ck.set(destination + cur, len);
			cur += len * sizeof(char);
			bm[tname] = ck;//store entries into the map
		}

		delete[] source;
		delete[] destination;
	}
	else //uncompressed data
	{
		fread(&NMap, sizeof(int32_t), 1, fp);
		for (int i = 0; i < NMap; ++i)
		{
			fread(cname, sizeof(char), 64, fp);
			string tname((char*)cname);
			fread(&ck.len, sizeof(int32_t), 1, fp);
			ck.p = new char[ck.len];
			if (NULL == ck.p) return false;
			fread(ck.p, sizeof(char), ck.len, fp);
			bm[tname] = ck;
		}
	}

	fclose(fp);
	return true;
}

double BinaryFile::getVersion()
{
	return header.version;
}
bool BinaryFile::write(const char* fname, bool b_compress, int32_t compressMethod)
{
	char cname[64];
	FILE* fp = fopen(fname, "wb");
	if (NULL == fp) return false;
	
	FileHeader header;
	header.MagicNum = BinaryFile_MagicNum; //indicate it's my file format
	header.b_compress = b_compress;
	header.compressMethod = compressMethod;
	header.version = BinaryFile_CurrentVerison;
	header.write(fp);
	int NMap = (int)bm.size();
	if (b_compress)
	{
		//count the data size first
		int32_t dataSize = sizeof(int32_t); //NMap
		for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
		{
			dataSize += sizeof(char) * 64 + sizeof(int32_t) + sizeof(char)*it->second.len;
		}
		char* source = new char[dataSize];
		char* destination = NULL;
		if (0 == compressMethod)
		{
			destination = new char[LZ4_compressBound(dataSize)];
		}
		else if (1 == compressMethod)
		{
			destination = new char[dataSize + 400];
		}
		//copy all the data to char* source
		size_t cur = 0;
		memcpy(source + cur, &NMap, sizeof(int));
		cur += sizeof(int);
		for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
		{
			if (it->first.size() > 63) return false;
			strcpy((char*)cname, it->first.c_str());
			memcpy(source + cur, cname, sizeof(char) * 64);
			cur += sizeof(char) * 64;
			memcpy(source + cur, &it->second.len, sizeof(int32_t));
			cur += sizeof(int32_t);
			memcpy(source + cur, it->second.p, sizeof(char)*it->second.len);
			cur += sizeof(char)*it->second.len;
		}
		//compress data from char* source to char* destination
		int compressed_size = 0;
		if (0 == compressMethod)
		{
			compressed_size = LZ4_compress(source, destination, dataSize);
		}
		else if (1 == compressMethod)
		{
			qlz_state_compress *qlz_state = new qlz_state_compress;
			int32_t compressed_size = (int32_t)qlz_compress(source, destination, dataSize, qlz_state);
			delete qlz_state;
		}
		//write the compressed data to file
		fwrite(&dataSize, sizeof(int32_t), 1, fp);
		fwrite(&compressed_size, sizeof(int32_t), 1, fp);
		fwrite(destination, compressed_size, 1, fp);
		delete[] source;
		delete[] destination;
	}
	else //don't need to compress the data
	{
		fwrite(&NMap, sizeof(int32_t), 1, fp);
		for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
		{
			if (it->first.size() > 63) return false;
			strcpy((char*)cname, it->first.c_str());
			fwrite(cname, sizeof(char), 64, fp);
			fwrite(&it->second.len, sizeof(int32_t), 1, fp);
			fwrite(it->second.p, sizeof(char), it->second.len, fp);
		}
	}
	fclose(fp);
	return true;
}
bool BinaryFile::write(FILE *fp) //write to steam without compression
{
	if (NULL == fp) return false;
	char* cname[64];
	int NMap = (int)bm.size();
	fwrite(&NMap, sizeof(int32_t), 1, fp);
	for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
	{
		if (it->first.size() > 63) return false;
		strcpy((char*)cname, it->first.c_str());
		fwrite(cname, sizeof(char), 64, fp);
		fwrite(&it->second.len, sizeof(int32_t), 1, fp);
		fwrite(it->second.p, sizeof(char), it->second.len, fp);
	}
	return true;
}
bool BinaryFile::writeChunk(const char* fname, string& key)
{
	if (bm.count(key) > 0)
	{
		char* cname[64];
		strcpy((char*)cname, key.c_str());

		void* p = bm[key].p;
		size_t len = bm[key].len;
		FILE* fp = fopen(fname, "wb");
		if (NULL == fp) return false;
		fwrite(cname, sizeof(char), 64, fp);
		fwrite(&len, sizeof(int32_t), 1, fp);
		fwrite(p, len, 1, fp);
		fclose(fp);
		return true;
	}
	else return false;
}
bool BinaryFile::get(const string& key, void* rp)
{
	if (bm.count(key) > 0)
	{
		Chunk& val = bm[key];
		memcpy(rp, val.p, val.len);
		return true;
	}
	return false; //no such a key
}
bool BinaryFile::getPointer(const string& key, void**p, int32_t &len)
{
	if (bm.count(key) > 0)
	{
		Chunk& val = bm[key];
		*p = val.p;
		len = val.len;
		return true;
	}
	*p = NULL;
	len = 0;
	return false; //no such a key
}
bool BinaryFile::push(const string& key, void* rp, int32_t len)
{
	Chunk ck;
	ck.len = len;
	ck.p = new char[len];
	memcpy(ck.p, rp, len);

	bm[key] = ck;
	return true;
}
bool BinaryFile::exist(const string& key)
{
	return bm.count(key) > 0;
}
int BinaryFile::size(){ return (int)bm.size(); }
void BinaryFile::beginSearch(){ ic = bm.begin(); }
bool BinaryFile::next(string& name, Chunk& chk)
{
	if (ic != bm.end())
	{
		name = ic->first;
		chk = ic->second;
		++ic;
		return true;
	}
	else return false;
}
bool BinaryFile::changeKeyName(string& oldKey, string& newKey)
{
	map<string, Chunk>::iterator it = bm.find(oldKey);
	if (it != bm.end())
	{
		if (bm.find(newKey) != bm.end()) return false; //the new key already exist
		bm[newKey] = it->second;
		bm.erase(oldKey);
		return true;
	}
	else return false;
}
bool BinaryFile::erase(const string& key)
{
	if (bm.find(key) != bm.end())
	{
		bm.erase(key);
		return true;
	}
	else return false;
}
bool BinaryFile::addChunk(string& key, Chunk& chk)
{
	if (bm.count(key) > 0) return false;
	bm[key] = chk;
	return true;
}
Chunk& BinaryFile::getChunk(const string& key)
{
	return bm[key];
}
int BinaryFile::getKeyLength(const string& key) //get key length in byte
{
	if (bm.count(key) > 0)
	{
		return bm[key].len;
	}
	return 0; //no such a key
}
int BinaryFile::getTotBytes()
{
	int sum = sizeof(int32_t);
	for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
	{
		sum += 64 + sizeof(int32_t);
		sum += it->second.len;
	}
	return sum;
}
void BinaryFile::delAll()
{
	for (map<string, Chunk>::iterator it = bm.begin(); it != bm.end(); ++it)
	{
		delete[] it->second.p;
	}
}
map<string, Chunk>& BinaryFile::getBM(){ return bm; }
/*******************************  End:  BinaryFile **************************************/

bool interp3N(ArrayMgr<SFloat>& in, ArrayMgr<SFloat>& out, int N)//interpolate 3D grid to N^3 more voxels
{
	if (1 == N)
	{
		out = in;
		return true;
	}
	int nx = in.getWidth(1);
	int ny = in.getWidth(2);
	int nz = in.getWidth(3);
	int MX = N * nx - N + 1;
	int MY = N * ny - N + 1;
	int MZ = N * nz - N + 1;
	
	if(!out.resize(MX, MY, MZ)) return false;//failed to resize memory because the matrix size is too large

	//interpolate
	SFloat sy0, sy1, sz0, sz1;
	SFloat FN = SFloat(N);
	for (int ix = 0; ix < nx - 1; ++ix)
		for (int iy = 0; iy < ny - 1; ++iy)
			for (int iz = 0; iz < nz - 1; ++iz)
			{
				for (int jx = 0; jx <= N; ++jx)
					for (int jy = 0; jy <= N; ++jy)
						for (int jz = 0; jz <= N; ++jz)
						{
							//interpolate in the z=0 plane
							//**interpolate in the y=0 line
							sy0 = in.a(ix, iy, iz)*(1 - jx / FN) + in.a(ix + 1, iy, iz)*(jx / FN);
							//**interpolate in the y=DY line
							sy1 = in.a(ix, iy + 1, iz)*(1 - jx / FN) + in.a(ix + 1, iy + 1, iz)*(jx / FN);
							sz0 = sy0*(1 - jy / FN) + sy1*(jy / FN);

							//interpolate in the z= DZ plane
							//**interpolate in the y=0 line
							sy0 = in.a(ix, iy, iz + 1)*(1 - jx / FN) + in.a(ix + 1, iy, iz + 1)*(jx / FN);
							//**interpolate in the y=DY line
							sy1 = in.a(ix, iy + 1, iz + 1)*(1 - jx / FN) + in.a(ix + 1, iy + 1, iz + 1)*(jx / FN);
							sz1 = sy0*(1 - jy / FN) + sy1*(jy / FN);

							//final interpolation along z direction
							out.a(N * ix + jx, N * iy + jy, N * iz + jz) = sz0*(1 - jz / FN) + sz1*(jz / FN);
						}
			}
	return true;
}

int gammaAnalysis(ArrayMgr<SFloat>& d1, ArrayMgr<SFloat>& d2, ArrayMgr<SFloat>& g, double dx, double dy, double dz, double percent, double DTA, bool local, int NInterp)
{
	//due to memory issue
	if (d1.getInnerLength() * sizeof(SFloat) * NInterp*NInterp*NInterp/ 1e9 > 1.5) return gammaAnalysis2(d1, d2, g, dx, dy, dz, percent, DTA, local, NInterp);

	ArrayMgr<SFloat> di1;
	if (NInterp > 1 && !interp3N(d1, di1, NInterp))	return -1;


	std::atomic<int> nPass(0);
	//DTA: distance to agree
	if (!d1.equalDims(d2)) return -1; //matrix dimension error

	int NSX = int(DTA / dx * NInterp);
	int NSY = int(DTA / dy * NInterp);
	int NSZ = int(DTA / dz * NInterp);

	if (0 == DTA) DTA = 0.01; //as denominator, it shouldn't be zero
	double DTA2 = DTA*DTA;
	int DNX = d1.getWidth(1);
	int DNY = d1.getWidth(2);
	int DNZ = d1.getWidth(3);

	int SK = NInterp;
	int SB = 1 - NInterp;
	//d2 is treated as basement, and scan d1 to find matches
	SFloat maxv, minv;
	d2.getMaxMin(maxv, minv);

#pragma omp parallel for
	for (int ix = 0; ix < DNX; ++ix) //iterate all the voxel of d2
		for (int iy = 0; iy < DNY; ++iy)
			for (int iz = 0; iz < DNZ; ++iz)
			{
				double deltaDose2 = (percent*maxv);
				deltaDose2 *= deltaDose2;
				SFloat min_g = 1e9;
				//search around this voxel from d1
				for (int jx = -NSX; jx <= NSX; ++jx)
					for (int jy = -NSY; jy <= NSY; ++jy)
						for (int jz = -NSZ; jz <= NSZ; ++jz)
						{
							//get the index of the searched voxel
							int six = SK*ix + jx;
							int siy = SK*iy + jy;
							int siz = SK*iz + jz;
							if (0 <= six&&six < SK*DNX + SB && 0 <= siy&&siy < SK*DNY + SB && 0 <= siz&&siz < SK*DNZ + SB)  //index inside the matrix
							{
								double dist2 = (jx*dx*jx*dx + jy*dy*jy*dy + jz*dz*jz*dz) / (SK*SK);
								if (dist2 < DTA2) //inside the sphere
								{
									double diff2 = NInterp > 1 ? di1.a(six, siy, siz) : d1.a(six, siy, siz);
									diff2 -= d2.a(ix, iy, iz);
									diff2 = diff2*diff2;

									if (local)
									{
										deltaDose2 = percent*d2.a(ix, iy, iz);
										deltaDose2 *= deltaDose2;
									}
									SFloat gi = (SFloat)sqrt(dist2 / DTA2 + diff2 / deltaDose2);
									if (gi < min_g) min_g = gi;
								}
							}
						}
				//found the min_g  => count pass number
				if (min_g <= 1)	++nPass;
				g.a(ix, iy, iz) = min_g;
			}

	return nPass;
}

int gammaAnalysis2(ArrayMgr<SFloat>& d1, ArrayMgr<SFloat>& d2, ArrayMgr<SFloat>& g, double dx, double dy, double dz, double percent, double DTA, bool local, int NInterp)
{
	//this function is designed for extreme large matrix, but has much slower speed than gammaAnalysis
	std::atomic<int> nPass(0);
	//DTA: distance to agree
	if (!d1.equalDims(d2)) return -1; //matrix dimension error

	int cx = (int)ceil(DTA / dx);
	int cy = (int)ceil(DTA / dy);
	int cz = (int)ceil(DTA / dz);


	int NSX = int(DTA / dx * NInterp);
	int NSY = int(DTA / dy * NInterp);
	int NSZ = int(DTA / dz * NInterp);

	if (0 == DTA) DTA = 0.01; //as denominator, it shouldn't be zero
	double DTA2 = DTA*DTA;
	int NDX = d2.getWidth(1);
	int NDY = d2.getWidth(2);
	int NDZ = d2.getWidth(3);

	//d2 is treated as basement, and scan d1 to find matches
	SFloat maxv, minv;
	d2.getMaxMin(maxv, minv);

#pragma omp parallel for
	for (int ix = 0; ix < NDX; ++ix) //iterate all the voxel of d2
		for (int iy = 0; iy < NDY; ++iy)
			for (int iz = 0; iz < NDZ; ++iz)
			{
				ArrayMgr<SFloat> d0(2 * cx + 1, 2 * cy + 1, 2 * cz + 1);
				ArrayMgr<SFloat> di(2 * NInterp*cx + 1, 2 * NInterp*cy + 1, 2 * NInterp*cy + 1);

				double deltaDose2 = (percent*maxv);
				deltaDose2 *= deltaDose2;
				SFloat min_g = 1e9;

				//need to prepare a group of voxels around (ix, iy, iz);
				//fetch the original voxels around(ix,iy,iz)
				for (int jx = -cx; jx <= cx; ++jx)
					for (int jy = -cy; jy <= cy; ++jy)
						for (int jz = -cz; jz <= cz; ++jz)
						{
							int ifx = ix + jx;
							int ify = iy + jy;
							int ifz = iz + jz;
							if (0 <= ifx&&ifx < NDX && 0 <= ify&&ify < NDY && 0 <= ifz&&ifz < NDZ)
							{
								d0.a(cx + jx, cy + jy, cz + jz) = d1.a(ifx, ify, ifz);
							}
							else d0.a(cx + jx, cy + jy, cz + jz) = 0;
						}
				interp3N(d0, di, NInterp);

				//search around this voxel from d1
				for (int jx = -NSX; jx <= NSX; ++jx)
					for (int jy = -NSY; jy <= NSY; ++jy)
						for (int jz = -NSZ; jz <= NSZ; ++jz)
						{
							//get the index of the searched voxel
							int six = NInterp*cx + jx;
							int siy = NInterp*cx + jy;
							int siz = NInterp*cx + jz;

							double dist2 = (jx*dx*jx*dx + jy*dy*jy*dy + jz*dz*jz*dz) / (NInterp*NInterp);
							if (dist2 < DTA2) //inside the sphere
							{
								double diff2 = di.a(six, siy, siz) - d2.a(ix, iy, iz);
								diff2 = diff2*diff2;
								if (local)
								{
									deltaDose2 = percent*d2.a(ix, iy, iz);
									deltaDose2 *= deltaDose2;
								}
								SFloat gi = (SFloat)sqrt(dist2 / DTA2 + diff2 / deltaDose2);
								if (gi < min_g) min_g = gi;
							}
						}
				//found the min_g  => count pass number
				if (min_g <= 1)	++nPass;
				g.a(ix, iy, iz) = min_g;
			}
	return nPass;
}
