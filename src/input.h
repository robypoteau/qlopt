#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include <misc.h>

namespace thesis{
	class input{
		int argc;
		char **argv;

		int extract(char* arg);
		vec stringToVector(string str);
		vec commaStringToVector(string str);

	public:
		input(){} //temp for testing
		input(int arg1, char** arg2){
			argc = arg1;
			argv = arg2;
		}
		string getSystem();
		vec getTimeData();
		vec getInterval();
		vec getU();
		vec getUGuess();
		vec getUNot();
		vec getYNot();
		bool isRegularized();
		bool isNoisy();
		double getNoise();
		int getNcoeffs();
		int getNumDivs();
		int getNumberOfIterations();
		bool useBSpline();
	};

	int input::extract(char* arg)
	{
		for(int i=0; i<argc; i++ ){
			if(strcmp(argv[i], arg) == 0){
				return i+1;
			}
		}
		return -1;
	}

	vec input::stringToVector(string str){

		if(/*type ==*/ true){
			size_t colon[3];
			colon[0] = 0;
			colon[1] = str.find(":");
			colon[2] = str.substr(colon[1]+1).find(":");

			double n[3];
			n[0]=strtod(str.substr(colon[0], colon[0+1]).c_str(), NULL);
			n[1]=strtod(str.substr(colon[1]+1, colon[1+1]).c_str(), NULL);
			n[2]=strtod(str.substr(colon[1] + colon[2] + 2).c_str(), NULL);

			int N = (int)((n[2]-n[0])/n[1] + 1.5);
			vec ans(N);
			ans(0) = n[0];
			for(int i=1; i<N; i++){
				ans(i) = ans(i-1) + n[1];
			}
			return ans;
		}
		else{
			vec ans;
			ans << 0;
			return ans;
		}
	}

	vec input::commaStringToVector(string str)
	{
		//...
		int counter = 0;
		size_t position = str.find(',');
		size_t changePosition = 0;
		std::vector<size_t> comma;
		comma.push_back(0);

		while(position != string::npos){
			counter++;
			changePosition = changePosition + position;
			position = str.substr(changePosition+counter).find(',');
			comma.push_back(changePosition+counter);
		}

		vec ans(counter+1);
		for(int i =0; i<counter; i++){
			ans(i) = strtod(str.substr(comma[i], comma[i+1]).c_str(), NULL);
		}
		ans(counter) = strtod(str.substr(comma[counter]).c_str(), NULL);

		return ans;
	}

	bool input::isRegularized()
	{
		if(extract((char*)"-r") != -1){
			return true;
		}
		else{
			return false;
		}
	}

	bool input::useBSpline(){
		if(extract((char*)"-b") != -1){
			return true;
		}
		else{
			return false;
		}
	}

	double input::getNoise()
	{
		return strtod(argv[extract((char*)"-n")],NULL);
	}

	int input::getNcoeffs()
	{
		return strtol(argv[extract((char*)"-p")],NULL,0);
	}

	int input::getNumDivs()
	{
		if(extract((char*)"-d") == -1){
			return 1;
		}else{
			return strtol(argv[extract((char*)"-d")],NULL,0);
		}
	}

	bool input::isNoisy(){
		if(extract((char*)"-n") != -1){
			return true;
		}
		else{
			return false;
		}
	}

	int input::getNumberOfIterations()
	{
		return strtol(argv[extract((char*)"-k")],NULL,0);
	}

	string input::getSystem()
	{
		return argv[extract((char*)"-s")];
	}

	vec input::getTimeData()
	{
		return stringToVector(argv[extract((char*)"-t")]);
	}

	vec input::getInterval()
	{
		return stringToVector(argv[extract((char*)"-i")]);
	}

	vec input::getU()
	{
		return commaStringToVector(argv[extract((char*)"-u")]);
	}

	vec input::getUGuess()
	{
		return commaStringToVector(argv[extract((char*)"-g")]);
	}

	vec input::getUNot()
	{
		return commaStringToVector(argv[extract((char*)"-o")]);
	}

	vec input::getYNot()
	{
		return commaStringToVector(argv[extract((char*)"-y")]);
	}
}
#endif
