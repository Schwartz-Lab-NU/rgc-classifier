#include <iostream>
#include "treeObjects.h"
#include <fstream>
#include <string>
#include <filesystem>

//usage: ./print PARAMS_ID FOLD# FOREST# TREE#
int main(int argc, char* argv[]) {

	std::filesystem::path treefile = (std::string)argv[1]
		+ "/fold" + (std::string) argv[2]
		+ "/forest" + (std::string) argv[3]
		+ "/tree" + (std::string) argv[4] +".out";

	std::cout << "Printing tree from ./" << treefile << std::endl;

	std::ifstream treestream(treefile, std::ios::binary);
	Tree tree(treestream);
	treestream.close();
	tree.print(false, false);

	return 0;
}