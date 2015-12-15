#include "C_tool.h"

C_tool::C_tool()
{

}

C_tool::~C_tool()
{
}

//void C_tool::debug(double tagNumber)
//{
//	std::cout<< tagNumber <<std::endl;
//	system("PAUSE");
//}


//void C_tool::getDir(std::string dirString, std::vector<std::string> & f)
//{
//	std::string newString = myreplace(dirString, "/", "\\");

//	FILE* pipe = NULL;
//	std::string pCmd = "dir /B /S " + newString;
//	char buf[256];

//	if (NULL == (pipe = _popen(pCmd.c_str(), "rt")))
//	{
//		cout << "Shit" << endl;
//		return;
//	}
//	while (!feof(pipe))
//	{
//		if (fgets(buf, 256, pipe) != NULL)
//		{
//			std::string str = std::string(buf);
//			str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
//			f.push_back(std::string(str));
//		}

//	}
//	_pclose(pipe);
//}

std::string C_tool::getFolder(const std::string& str)
{
	size_t found,end;
	found = str.find_last_of("/\\");
	end = str.length();
	//cout << " folder: " << str.substr(0, found) << endl;
	//cout << " file: " << str.substr(found + 1) << endl;
	return str.substr(0, found);
}

std::string C_tool::getFileName(const std::string& str)
{
	size_t found;
	found = str.find_last_of("/\\");
	return str.substr(found + 1);
}

std::string C_tool::getExtension(std::string str)
{
	size_t found;
	found = str.find_last_of(".");
	std::string returnString = str.substr(found + 1, 3);
	return returnString;
}

bool C_tool::isDicom(std::string str)
{
	std::string extension = getExtension(str);
	return ((extension == "dcm") ? true : false) ;
}

//__int64 C_tool::GetFileSize(std::wstring const &path) {

//	WIN32_FIND_DATAW data;
//	HANDLE h = FindFirstFileW(path.c_str(), &data);
//	if (h == INVALID_HANDLE_VALUE)
//		return -1;

//	FindClose(h);

//	return data.nFileSizeLow | (__int64)data.nFileSizeHigh << 32;
//}

std::string C_tool::myreplace(std::string &subject, const std::string &search, const std::string &replace)
{
	size_t pos = 0;
	while ((pos = subject.find(search, pos)) != std::string::npos) {
		subject.replace(pos, search.length(), replace);
		pos += replace.length();
	}
	return subject;

	//return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

//void C_tool::makeFolder(std::string aPath)
//{
//	char DirName[256];
//	const char* p = aPath.c_str();
//	char* q = DirName;
//	while (*p)
//	{
//		if (('\\' == *p) || ('/' == *p))
//		{
//			if (':' != *(p - 1))
//			{
//				CreateDirectory(DirName, NULL);
//			}
//		}
//		*q++ = *p++;
//		*q = '\0';
//	}
//	CreateDirectory(DirName, NULL);
//}



void C_tool::showSize(signedShort3D::Pointer inputImage)
{
	signedShort3D::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
	signedShort3D::SizeType aSize = inputRegion.GetSize();
	cout << aSize[0] << "," << aSize[1] << "," << aSize[2] << endl;
}

void C_tool::showOrigin(signedShort3D::Pointer inputImage)
{
	signedShort3D::PointType origin = inputImage->GetOrigin();
	cout << origin[0] << "," << origin[1] << "," << origin[2] << endl;
}

void C_tool::showIndexStart(signedShort3D::Pointer inputImage)
{
	signedShort3D::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
	signedShort3D::IndexType start = inputRegion.GetIndex();
	cout << start[0] << "," << start[1] << "," << start[2] << endl;
}


