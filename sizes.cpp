#include <iostream>
#include "treeObjects.h"

int main(int argc, char* argv[]) {
  std::cout << "Sizeof:" << std::endl;
  std::cout << "\tParamSet: " << sizeof(paramSet) << std::endl;
  std::cout << "\tRandomDraws: " << sizeof(randomDraws) << std::endl;
	return 0;
}