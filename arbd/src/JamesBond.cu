/* Bond, James Bond.cu Copyright Justin Dufresne and Terrance Howard, 2013.
 * We prefer our code shaken, not stirred.
 */

#include "JamesBond.h"

Bond::Bond(String strflag, int ind1, int ind2, String fileName) :
		ind1(ind1), ind2(ind2), fileName(fileName) {
	if (strflag == "REPLACE") {
		flag = REPLACE;
	} else if (strflag == "ADD") {
		flag = ADD;
	} else {
		printf("WARNING: Invalid operation flag found:"
					 "\"BOND %s %d %d\"\n", strflag.val(), ind1, ind2);
		printf("         Using default flag\n");
		flag = DEFAULT;
	}
	tabFileIndex = -1;
}

void Bond::print() {
	printf("BOND %s %d %d %s\n", flags[flag].val(), ind1, ind2, fileName.val());
}

String Bond::toString() {
	return "BOND " + flags[flag] + " " + ind1 + " " + ind2 + " " + fileName;
}
