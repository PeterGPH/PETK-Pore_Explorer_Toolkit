#ifndef USEFUL_CU
#define USEFUL_CU

///////////////////////////////////////////////////////////////////////
// Useful classes.
// Author: Jeff Comer <jcomer2@illinois.edu>
#include "useful.h"
#include <cuda.h>


bool isReal(char c) {
	char num[] = "0123456789-eE.";

	for (int i = 0; i < 14; i++) {
		if (c == num[i]) return true;
	}
	return false;
}

bool isInt(char c) {
	char num[] = "0123456789-";

	for (int i = 0; i < 11; i++) {
		if (c == num[i]) return true;
	}
	return false;
}

int firstSpace(const char* s, int max) {
	for (int i = 0; i < max; i++) {
		if (s[i] == ' ') return i;
	}
	return -1;
}

// A classic growable string class.

void String::print() const {
	printInline();
	printf("\n");
}
void String::printInline() const {
	for (int i = 0; i < len; i++)
		printf("%c", c[i]);
}

String& String::operator=(const String& s) {
	len = s.len;
	cap = s.len;
	c = new char[cap];
	for (int i = 0; i < s.len; i++) {
		c[i] = s.c[i];
	}
	return *this;
}

void String::add(char c0) {
	int n = 2;
	len--;
	if (n + len > cap) grow(n + len);
	c[len] = c0;
	c[len+1] = '\0';
	len += n;
}

void String::add(const char* s) {
	const int n = strlen(s) + 1;
	len--;
	if (n + len > cap) grow(n + len);
	for (int i = 0; i < n; i++) c[i+len] = s[i];
	len += n;
}

void String::add(String& s) {
	len--;

	if (len + s.len > cap) grow(len + s.len);
	for (int i = 0; i < s.len; i++) c[i+len] = s.c[i];
	len += s.len;
}

int String::length() const {
	return len-1;
}

// Negative indices go from the end.
String String::range(int first, int last) const {
	String ret;
	const int len1 = len - 1;

	first = first % len1;
	last = last % len1;

	if (first < 0) first += len1;
	if (last < 0) last += len1;
	for(int i = first; i <= last; i++) ret.add(c[i]);

	return ret;
}

String String::trim() {
	int i0, i1;
	for (i0 = 0; i0 < len; i0++) if (!isWhite(c[i0])) break;
	for (i1 = len-1; i1 >= 0; i1++) if (!isWhite(c[i0])) break;
	return range(i0,i1);
}

bool String::isWhite(char c) {
	char num[] = " \n\t\v\b\r\f\a";

	for (int i = 0; i < 8; i++) {
		if (c == num[i]) return true;
	}
	return false;
}

int String::tokenCount() const {
	const int len1 = len - 1;

	int count = 0;

	// Skip leading whitespace.
	int i = 0;
	while (isWhite(c[i]) && i < len1) i++;

	if (i >= len1) return 0;

	// Find the tokens.
	bool onToken = true;
	count++;
	while (i < len1) {
		if (onToken) {
			if (isWhite(c[i])) onToken = false;
		} else {
			if (!isWhite(c[i])) {
				onToken = true;
				count++;
			}
		}
		i++;
	}

	return count;
}

int String::tokenize(String* tokenList) const {
	const int len1 = len - 1;
	const int n = tokenCount();

	// Skip leading whitespace.
	int i = 0;
	while (isWhite(c[i]) && i < len1) i++;

	// Find the tokens.
	bool onToken = true;
	int j = 0;
	int tokenStart = i;
	while (i < len1) {
		if (onToken) {
			if (isWhite(c[i])) {
				onToken = false;
				tokenList[j] = range(tokenStart,i-1);
				//printf("(%s)\n", tokenList[j].val());
				j++;
			}
		} else if (!isWhite(c[i])) {
			onToken = true;
			tokenStart = i;
		}
		i++;
	}

	// Grab the last token if there is one.
	if (onToken) tokenList[n-1] = range(tokenStart,-1);
	return n;
}

int String::tokenCount(char delim) const {
	const int len1 = len - 1;

	int count = 0;

	// Skip leading whitespace.
	int i = 0;
	while (c[i] == delim && i < len1) i++;

	if (i >= len1) return 0;

	// Find the tokens.
	bool onToken = true;
	count++;
	while (i < len1) {
		if (onToken) {
			if (c[i] == delim) onToken = false;
		} else if (c[i] != delim) {
			onToken = true;
			count++;
		}
		i++;
	}

	return count;
}

int String::tokenize(String* tokenList, char delim) const {
	const int len1 = len - 1;
	const int n = tokenCount(delim);

	// Skip leading whitespace.
	int i = 0;
	while (c[i] == delim && i < len1) i++;

	// Find the tokens.
	bool onToken = true;
	int j = 0;
	int tokenStart = i;
	while (i < len1) {
		if (onToken) {
			if (c[i] == delim) {
				onToken = false;
				tokenList[j] = range(tokenStart,i-1);
				//printf("(%s)\n", tokenList[j].val());
				j++;
			}
		} else {
			if (c[i] != delim) {
				onToken = true;
				tokenStart = i;
			}
		}
		i++;
	}

	// Grab the last token if there is one.
	if (onToken) tokenList[n-1] = range(tokenStart,-1);
	return n;
}

void String::lower() {
	for(int i = 0; i < len; i++) {
		if (c[i] >= 0x41 && c[i] <= 0x5A)
			c[i] += 0x20;
	}
}

void String::upper() {
	for(int i = 0; i < len; i++) {
		if (c[i] >= 0x61 && c[i] <= 0x7A)
			c[i] -= 0x20;
	}
}

char String::operator[](int i) const {
	if (i < 0 || i >= len) return '\0';
	return c[i];
}

bool String::operator==(const String& s) const {
	if (len != s.len) return false;
	for (int i = 0; i < len; i++)
		if (c[i] != s.c[i])
			return false;
	return true;
}

bool String::operator==(const char* s) const {
	return (*this == String(s));
}

bool String::operator!=(const String& s) const {
	return !(*this == s);
}

const char* String::val() const {
	return c;
}

String String::getNumbers() const {
	String ret;
	for(int i = 0; i < len; i++) {
		if (isInt(c[i])) ret.add(c[i]);
	}
	return ret;
}

void String::grow(int n) {
	char* c0 = c;
	c = new char[n];
	cap = n;
	for (int i = 0; i < len; i++) c[i] = c0[i];
	delete[] c0;
}

String operator+(String s, int i) {
	String ret("");
	ret.add(s);
	char* c = new char[256];
	int digit = 0;
	while (i > 0)
	{
		char next = char(i % 10 + (int)'0');
		c[digit++] = next;
		i /= 10;
	}
	for (int j = digit - 1; j >= 0; j--)
	{
		ret.add(c[j]);
	}
	delete[] c;
	return ret;
}

String operator+(String s, const char* c) {
	String ret("");
	ret.add(s);
	ret.add(c);	
	return ret;
}

String operator+(String s1, String s2) {
	String ret("");
	ret.add(s1);
	ret.add(s2);
	return ret;
}

// class Vector3
// Operations on 3D float vectors
//
Vector3 Vector3::random(float s) {
	Vector3 v;
	v.x = (float(rand())/RAND_MAX-0.5f)*s;
	v.y = (float(rand())/RAND_MAX-0.5f)*s;
	v.z = (float(rand())/RAND_MAX-0.5f)*s;
	return v;
}

String Vector3::toString() const {
	char s[128];
	sprintf(s, "%.10g %.10g %.10g", x, y, z);
	return String(s);
}

// class Matrix3
// Operations on 3D float matrices
//
HOST DEVICE	
Matrix3::Matrix3(float s) {
	exx = s;
	exy = 0.0f;
	exz = 0.0f;
	eyx = 0.0f;
	eyy = s;
	eyz = 0.0f;
	ezx = 0.0f;
	ezy = 0.0f;
	ezz = s;
	isDiag = true;
}

Matrix3::Matrix3(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz) {
	exx = xx;
	exy = xy;
	exz = xz;
	eyx = yx;
	eyy = yy;
	eyz = yz;
	ezx = zx;
	ezy = zy;
	ezz = zz;
	setIsDiag();
}

Matrix3::Matrix3(float x, float y, float z) {
	exx = x;
	exy = 0.0f;
	exz = 0.0f;
	eyx = 0.0f;
	eyy = y;
	eyz = 0.0f;
	ezx = 0.0f;
	ezy = 0.0f;
	ezz = z;
	isDiag = true;
}

Matrix3::Matrix3(const Vector3& ex, const Vector3& ey, const Vector3& ez) {
	exx = ex.x;
	eyx = ex.y;
	ezx = ex.z;
	exy = ey.x;
	eyy = ey.y;
	ezy = ey.z;
	exz = ez.x;
	eyz = ez.y;
	ezz = ez.z;
	setIsDiag();
}

Matrix3::Matrix3(const float* d) {
	exx = d[0];
	exy = d[1];
	exz = d[2];
	eyx = d[3];
	eyy = d[4];
	eyz = d[5];
	ezx = d[6];
	ezy = d[7];
	ezz = d[8];	
	setIsDiag();
}	
Matrix3 Matrix3::inverse() const {
	Matrix3 m;
	if (isDiag) {
		m = Matrix3(1.0f/exx,1.0f/eyy,1.0f/ezz);
	} else {
		float det = exx*(eyy*ezz-eyz*ezy) - exy*(eyx*ezz-eyz*ezx) + exz*(eyx*ezy-eyy*ezx);
		m.exx = (eyy*ezz - eyz*ezy)/det;
		m.exy = -(exy*ezz - exz*ezy)/det;
		m.exz = (exy*eyz - exz*eyy)/det;
		m.eyx = -(eyx*ezz - eyz*ezx)/det;
		m.eyy = (exx*ezz - exz*ezx)/det;
		m.eyz = -(exx*eyz - exz*eyx)/det;
		m.ezx = (eyx*ezy - eyy*ezx)/det;
		m.ezy = -(exx*ezy - exy*ezx)/det;
		m.ezz = (exx*eyy - exy*eyx)/det;
		m.isDiag = isDiag;
	}
	return m;
}

float Matrix3::det() const {
	return isDiag ? exx*eyy*ezz :
		exx*(eyy*ezz-eyz*ezy) - exy*(eyx*ezz-eyz*ezx) + exz*(eyx*ezy-eyy*ezx);
}

String Matrix3::toString() const {
	char s[128];
	sprintf(s, "%2.8f %2.8f %2.8f\n%2.8f %2.8f %2.8f\n%2.8f %2.8f %2.8f",
					exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz);
	return String(s);
}

String Matrix3::toString1() const {
	char s[128];
	sprintf(s, "%2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f %2.8f",
					exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz);
	return String(s);
}


// class IndexList
// A growable list of integers.
IndexList::IndexList() {
	num = 0;
	maxnum = 16;
	lis = new int[maxnum];
}

IndexList::IndexList(const IndexList& l) {
	num = l.num;
	maxnum = num + 16;
	lis = new int[maxnum];

	for(int i = 0; i < l.num; i++) lis[i] = l.lis[i];
}

IndexList::IndexList(const int* a, int n) {
	num = n;
	maxnum = num + 16;
	lis = new int[maxnum];

	for(int i = 0; i < n; i++) lis[i] = a[i];
}

IndexList::~IndexList() {
	delete[] lis;
}

void IndexList::add(const int val) {
	// If we need more space, allocate a new block that is 1.5 times larger
	// and copy everything over
	if (num == maxnum) {
		maxnum = (maxnum*3)/2 + 1;
		int* oldlis = lis;
		lis = new int[maxnum];
		int i;
		for(i = 0; i < num; i++)
			lis[i] = oldlis[i];
		
		delete [] oldlis;
	}

	// We should have enough space now, add the value
	lis[num] = val;
	num++;
}

IndexList& IndexList::operator=(const IndexList& l) {
	delete[] lis;

	num = l.num;
	maxnum = num + 16;
	lis = new int[maxnum];

	for(int i = 0; i < num; i++) lis[i] = l.lis[i];
	return *this;
}

void IndexList::add(const IndexList& l) {
	int oldnum = num;
	num = num + l.num;

	if (num > maxnum) {
		maxnum = (num*3)/2 + 1;
		int* oldlis = lis;
		lis = new int[maxnum];

		for(int i = 0; i < oldnum; i++) lis[i] = oldlis[i];
		delete[] oldlis;
	}

	for(int i = 0; i < l.num; i++) lis[i+oldnum] = l.lis[i];
}

int* IndexList::getList() {
	return lis;
}

void IndexList::clear() {
	num=0;
	maxnum=16;
	delete[] lis;
	lis = new int[maxnum];
}

String IndexList::toString() const {
	String ret;
	char tmp[32];

	for (int i = 0; i < num; i++) {
		sprintf(tmp, "%i ", lis[i]);
		ret.add(tmp);
	}
	return ret;
}

int IndexList::find(int key) {
	for(int i = 0; i < num; i++)
		if (lis[i] == key)
			return i;
	return -1;
}

IndexList IndexList::range(int i0, int i1) {
	if (i0 >= i1)
		return IndexList();

	if (i0 < 0) i0 = 0;
	if (i1 >= num) i1 = num-1;
	int n = i1 - i0 + 1;

	IndexList ret(&(lis[i0]), n);
	return ret;
}

#endif
