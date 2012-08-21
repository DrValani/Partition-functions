#ifndef Storage_CC
#define Storage_CC
#include <vector>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include "latticeSize.cc"
#include "Globals.cc"
#include "Polynomial.cc"
using namespace std;
// {{{ Storage class

//class Storage;
//ostream & operator<<(ostream & out, const Storage &s);

/**
 *   For storing vlists the size is 2^xsize
 */
class Storage {
public:

	Storage() :
			group_id(-1) {
	}
	Storage(int a);  //:group_id(-1),group_id2(-1);

	Storage(bool _clear) :
			group_id(-1) {
		if (_clear) clear();
	}  // Create storage with -1 entries
	void push_back(const int index, int pow);

	int getPow(const int i) const {
		return list[i];
	}

	void setPow(const int index, const int i) {
		list[index] = i;
	}
	void print() const;
	void clear();

	bool operator==(const Storage &s) const;
	bool operator<(const Storage &s) const;
	bool equals(const Storage &s) const;  // compare everything

	int operator[](const int i) const {
		return list[i];
	}
	Storage order_it() const;
	Storage order_it_asending() const;
	void reverse(Storage temp);

	int getId() const {
		return group_id;
	}

	void setId(int id) {
		group_id = id;
	}

	void sum(const vector<Polynomial> &vpol, Polynomial &p) const;
	template<class T>
	T applyGroup(const T &s1) const;
	bool check_structure() const;
	void operator+(const int &x);

	void createGroup(const Storage &s1, const Storage &s2);

	int searchGroup(vector<Storage> &allg1) const;
private:
	int list[ARRSIZE];
	int group_id;
};

void Storage::createGroup(const Storage &s1, const Storage &s2) {
	for (int i = 0; i < ARRSIZE; i++) {
		int j = ARRSIZE - 1;
		while (j >= 0) {
			if (s1[j] == s2[i]) {
				list[i] = j;
				break;
			}
			j--;
		}
	}
}

int Storage::searchGroup(vector<Storage> &allg1) const {

	unsigned int i = 0;
	while (i < allg1.size()) {
		if (allg1[i].equals(*this)) {
			return i;
		}
		i++;
	}
	allg1.push_back(*this);
	return allg1.size() - 1;
}

Storage::Storage(int a) :
		group_id(-1) {
	for (int i = 0; i < ARRSIZE; i++)
		list[i] = a;
}

void Storage::operator+(const int &x) {
	for (int i = 0; i < ARRSIZE; i++)
		list[i] = list[i] + x;
}

/**
 * This will check that every g[k] has at least one element
 * where k is 0 to ARRSIZE-1
 * and that every element(0 to ARRSIZE-1)  is present.
 **/
bool Storage::check_structure() const {
	Storage counter;
	counter.clear();
	for (int i = 0; i < ARRSIZE; i++) {
		counter.push_back(list[i], 1);
	}
	for (int i = 0; i < ARRSIZE; i++)
		if (counter[i] == -1) return false;
	return true;
}

template<class T>
T Storage::applyGroup(const T &s1) const {
	T temp;
	for (int i = 0; i < ARRSIZE; i++) {
		int pow = s1[this->list[i]];
		temp.push_back(i, pow);
	}
	return temp;
}

void Storage::sum(const vector<Polynomial> &vpol, Polynomial &p) const {
	p.clear();
	int ct = 1;
	for (int i = 0; i < ARRSIZE - 1; i++) {
		if (list[i] == list[i + 1]) {
			ct++;
		} else {
			p.multiPlus(vpol[list[i]], ct);
			ct = 1;
		}
	}
	p.multiPlus(vpol[list[ARRSIZE - 1]], ct);
}

void Storage::reverse(Storage temp) {
	for (int i = 0; i < ARRSIZE; i++) {
		this->push_back(i, temp[ARRSIZE - 1 - i]);
	}
}

Storage Storage::order_it_asending() const {
	Storage temp;
	multiset<int> s(list, list + ARRSIZE);
	multiset<int>::iterator sittr, eittr;
	sittr = s.begin();
	eittr = s.end();
	int i = 0;
	for (; sittr != eittr; sittr++) {
		temp.push_back(i, *sittr);
		i++;
	}
	return temp;
}

Storage Storage::order_it() const {
	Storage temp;
	multiset<int> s(list, list + ARRSIZE);
	multiset<int>::iterator sittr, eittr;
	sittr = s.begin();
	eittr = s.end();
	int i = ARRSIZE;
	for (; sittr != eittr; sittr++) {
		i--;
		temp.push_back(i, *sittr);
	}
	return temp;
}

bool Storage::operator<(const Storage &s) const {
	for (int i = 0; i < ARRSIZE; i++) {
		if (this->list[i] < s[i]) return true;
		if (this->list[i] > s[i]) return false;
	}
	return false;
}

bool Storage::operator==(const Storage &s) const {
	Storage temp = this->order_it();
	Storage temp2 = s.order_it();
	return temp.equals(temp2);
}

bool Storage::equals(const Storage &s) const {
	for (int i = 0; i < ARRSIZE; i++)
		if (list[i] != s[i]) return false;
	return true;
}

void Storage::clear() {
	for (int i = 0; i < ARRSIZE; i++)
		list[i] = -1;
}

void Storage::push_back(const int index, int pow) {
	this->list[index] = pow;
}

void Storage::print() const {
	for (int i = 0; i < ARRSIZE; i++)
		cout << list[i] << " ";
	cout << endl;
}

// }}}

ostream & operator<<(ostream & out, const Storage &s) {
    for (int i = 0; i < ARRSIZE; i++) {
        out << s[i] << " ";
    }
    out << " -- groupid =" << s.getId();
    return out;
}

vector<Storage> TRANSFORMGROUP;
int search_add_list_binary(vector<Storage> &vs2, Storage &s2) {
    if (vs2.size() == 0) {
        vs2.push_back(s2);
        return 0;
    }
    int low = 0;
    int high = vs2.size() - 1;
    int middle = 0;
    while (low <= high) {
        middle = (low + high) / 2;
        if (vs2[middle].equals(s2))
            return middle;
        else if (vs2[middle] < s2) {
            high = middle - 1;
        } else {
            low = middle + 1;
        }
    }
    if (vs2[middle] < s2) {
        vs2.insert(vs2.begin() + middle, s2);
        return middle;
    }
    vs2.insert(vs2.begin() + middle + 1, s2);
    return middle + 1;

}


#endif

