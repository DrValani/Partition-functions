#ifndef Group_CC
#define Group_CC

#include <algorithm>
#include <set>
#include "Storage.cc"
#include "Globals.cc"
using namespace std;
// {{{ Group Class
class Group;
typedef vector<Group> vGroup;
typedef vector<vGroup> vvGroup;

/**
 * Groups Class
 */
class Group {
public:

	Group() :
			smallest_form(true) {
	}
	template<class T>
	void createGroup(T &s);
	template<class T>
	void createGroup(const T &s1, const T &s2);
	template<class T>
	bool createGroup1(const T &s1, const T &s2, const int x);
	void createSubGroups(vGroup &vg);
	void createSubGroups(vGroup &vg, int i);
	int searchGroup(vGroup &allg1) const;
	int searchGroup(vector<Storage> &subg) const;
	template<class T>
	int searchGroup(vGroup &allg1, T &subGroup) const;
	template<class T>
	T applyGroup(const T &s1) const;
	Storage applyInverse(const Storage &s1);
	void print() const;
	void print(const Storage &s) const;
	bool isSubset(Group &allg1, Group g1) const;
	Group createSubset(const Group &g1) const;
	bool operator==(const Group &s) const;

	int getOrder(int &i) const {
		return *g[i].begin();
	}

	int operator[](const int i) const {
		return *g[i].begin();
	}
	void push_back(int i, int pow);
	void remove(int j, int &k);
	bool check_stability() const;
	bool check_structure() const;
	bool check_smallestform() const;
	bool check_smallestform2();
	bool equals(const Group & gp) const;

	bool isSmallestform() const {
		return smallest_form;
	}
	void reduce();
	void convert(Storage &s) const;
private:
	set<int> g[ARRSIZE];
	bool smallest_form;  // Cannot reduce the group anymore
};

void Group::convert(Storage &s) const {
	for (int i = 0; i < ARRSIZE; i++) {
		s.push_back(i, *g[i].begin());
	}
}

void Group::reduce() {
	for (int i = 0; i < ARRSIZE; i++) {
		if (g[i].size() > 1) {
			set<int>::iterator s1 = g[i].begin();
			this->remove(*s1, i);
		}
	}
	if (!this->check_structure()) {
		cout << "Error in Group::reduce" << endl;
		exit(1);
	}

}

bool Group::equals(const Group & gp) const {
	for (int i = 0; i < ARRSIZE; i++) {
		if (this->g[i] != gp.g[i]) return false;
	}
	return true;
}

void Group::createSubGroups(vGroup &vg) {
	if (this->smallest_form) {
		this->searchGroup(vg);
	} else {
		for (int i = 0; i < ARRSIZE; i++) {
			if (g[i].size() > 1) {
				set<int>::iterator s1, e1;
				s1 = g[i].begin();
				e1 = g[i].end();
				for (; s1 != e1; s1++) {
					Group temp = *this;
					temp.remove(*s1, i);
					if (temp.check_structure()) {
						temp.smallest_form = temp.check_smallestform();
						temp.createSubGroups(vg, i);
					}
				}
			}
			cout << i << endl;
		}
	}
}

void Group::createSubGroups(vGroup &vg, int i) {
	if (this->smallest_form) {
		this->searchGroup(vg);
	} else {
		for (; i < ARRSIZE; i++) {
			if (g[i].size() > 1) {
				set<int>::iterator s1, e1;
				s1 = g[i].begin();
				e1 = g[i].end();
				for (; s1 != e1; s1++) {
					Group temp = *this;
					temp.remove(*s1, i);
					if (temp.check_structure()) {
						temp.smallest_form = temp.check_smallestform();
						temp.createSubGroups(vg, i);
					}
				}
			}
		}
	}
}

bool Group::check_smallestform() const {
	for (int k = 0; k < ARRSIZE; k++) {
		if (g[k].size() > 1) return false;
	}
	return true;
}

bool Group::check_smallestform2() {
	for (int k = 0; k < ARRSIZE; k++) {
		if (g[k].size() > 1) return false;
	}
	smallest_form = true;
	return true;
}

bool Group::check_stability() const {
	for (int k = 0; k < ARRSIZE; k++) {
		if (g[k].size() < 1) return false;
	}
	return true;
}

/**
 * This will check that every g[k] has at least one element
 * where k is 0 to ARRSIZE-1
 * and that every element(0 to ARRSIZE-1)  is present.
 **/
bool Group::check_structure() const {
	for (int k = 0; k < ARRSIZE; k++) {
		if (g[k].size() < 1) return false;
	}
	set<int>::iterator s1, e1;
	Storage counter;
	for (int i = 0; i < ARRSIZE; i++)
		counter.push_back(i, 0);
	for (int i = 0; i < ARRSIZE; i++) {
		s1 = g[i].begin();
		e1 = g[i].end();
		while (s1 != e1) {
			counter.push_back(*s1, 1);
			s1++;
		}
	}
	for (int i = 0; i < ARRSIZE; i++)
		if (counter[i] == 0) return false;

	return true;
}

/**
 * Remove all j, except the j in position i.
 */
void Group::remove(int j, int &i) {
	for (int k = 0; k < ARRSIZE; k++) {
		if (k != i) {
			set<int>::iterator st = this->g[k].begin();
			set<int>::iterator end = this->g[k].end();
			while (st != end) {
				if (*st == j) {
					this->g[k].erase(st);
					st = this->g[k].begin();
					end = this->g[k].end();
				} else st++;
			}
		} else {
			this->g[k].clear();
			this->g[k].insert(j);
		}
	}
}

Storage Group::applyInverse(const Storage &s1) {
	Storage temp;
	Group tempGroup = *this;
	if (!smallest_form) {
		// Then make into the smallest form
		for (int i = 0; i < ARRSIZE; i++) {
			int j = tempGroup.getOrder(i);
			temp.push_back(j, s1[i]);
			tempGroup.remove(j, i);
		}
		if (tempGroup.check_structure()) {
			Group g1;
			g1.createGroup(temp, s1);
			*this = this->createSubset(g1);
			return temp;
		} else {
			cout << "Cannot create inverse\n";
			tempGroup.print();
			this->print();
			s1.print();
			temp.print();
			exit(0);
		}
	}
	for (int i = 0; i < ARRSIZE; i++) {
		int j = this->getOrder(i);
		temp.push_back(j, s1[i]);
	}
	return temp;
}

void Group::push_back(int i, int pow) {
	g[i].insert(pow);
	if (g[i].size() > 1) smallest_form = false;
}

template<class T>
T Group::applyGroup(const T &s1) const {
	T temp;
	for (int i = 0; i < ARRSIZE; i++) {
		int pow = s1[this->getOrder(i)];
		temp.push_back(i, pow);
	}
	return temp;
}

int Group::searchGroup(vGroup &allg1) const {
	bool found = false;
	unsigned int i = 0;
	while (i < allg1.size() && !found) {
		found = allg1[i].isSubset(allg1[i], *this);
		if (found) {
			return i;
		}
		i++;
	}
	allg1.push_back(*this);
	return allg1.size() - 1;
}

int Group::searchGroup(vector<Storage> &allg1) const {
	Storage s;
	if (this->smallest_form) this->convert(s);
	else {
		Group g1 = *this;
		g1.reduce();
		g1.convert(s);
	}
	unsigned int i = 0;
	while (i < allg1.size()) {
		if (allg1[i].equals(s)) return i;
		i++;
	}
	allg1.push_back(s);
	return allg1.size() - 1;
}

template<class T>
int Group::searchGroup(vGroup &allg1, T &subGroup) const {
	int found = this->searchGroup(allg1);
	return allg1[found].searchGroup(subGroup);
}

Group Group::createSubset(const Group &g1) const {
	Group temp;
	set<int>::iterator sthis, s1, ethis, e1;
	for (int k = 0; k < ARRSIZE; k++) {
		sthis = g[k].begin();
		ethis = g[k].end();
		while (sthis != ethis) {
			s1 = g1.g[k].begin();
			e1 = g1.g[k].end();
			bool found1 = false;
			while (s1 != e1 && !found1) {
				if (*sthis == *s1) {
					temp.g[k].insert(*sthis);
					found1 = true;
				}
				s1++;
			}
			if (temp.g[k].size() > 1) temp.smallest_form = false;
			sthis++;
		}
	}
	return temp;
}

bool Group::operator==(const Group &s) const {
	for (int i = 0; i < ARRSIZE; i++)
		if (g[i] != s.g[i]) return false;
	return true;
}

template<class T>
void Group::createGroup(const T &s1, const T &s2) {
	bool found = false;
	for (int i = 0; i < ARRSIZE; i++) {
		int j = ARRSIZE - 1;
		while (j >= 0) {
			if (s1[j] == s2[i]) {
				g[i].insert(j);
				if (!found) if (g[i].size() > 1) {
					smallest_form = false;
					found = true;
				}
			}
			j--;
		}
	}
}

template<class T>
void Group::createGroup(T &s) {
	bool found = false;
	for (int i = 0; i < ARRSIZE; i++) {
		int j = 0;
		while (j < ARRSIZE) {
			if (s[j] == s[i]) {
				g[i].insert(j);
				if (!found) if (g[i].size() > 1) {
					smallest_form = false;
					found = true;
				}
			}
			j++;
		}
	}
}

bool Group::isSubset(Group &allg1, Group g1) const {
	Group temp;
	set<int>::iterator sthis, s1, ethis, e1;
	for (int k = 0; k < ARRSIZE; k++) {
		sthis = g[k].begin();
		ethis = g[k].end();
		while (sthis != ethis) {
			s1 = g1.g[k].begin();
			e1 = g1.g[k].end();
			while (s1 != e1) {
				if (*sthis == *s1) {
					temp.g[k].insert(*s1);
					s1 = e1;
				} else s1++;
			}
			sthis++;
		}
		if (temp.g[k].empty()) return false;
		else if (temp.g[k].size() == 1) {
			g1.remove(temp.getOrder(k), k);
		} else if (temp.g[k].size() > 1) temp.smallest_form = false;
	}
	// Make sure every element is still there
	Storage counter;
	for (int i = 0; i < ARRSIZE; i++)
		counter.push_back(i, 0);
	for (int i = 0; i < ARRSIZE; i++) {
		s1 = temp.g[i].begin();
		e1 = temp.g[i].end();
		while (s1 != e1) {
			counter.push_back(*s1, 1);
			s1++;
		}
	}
	for (int i = 0; i < ARRSIZE; i++)
		if (counter[i] == 0) return false;
	allg1 = temp;
	return true;
}

void Group::print() const {
	for (int k = 0; k < ARRSIZE; k++) {
		set<int>::iterator sitrr, eitrr;
		sitrr = g[k].begin();
		eitrr = g[k].end();
		// cout << "(";
		while (sitrr != eitrr) {
			cout << *sitrr << " ";
			sitrr++;
		}
		cout << ", ";
	}
	if (smallest_form) cout << "irreduceble";
	else cout << "can be reduced";
	cout << endl;
}

void Group::print(const Storage &s) const {
	for (int k = 0; k < ARRSIZE; k++) {
		cout << s[getOrder(k)] << " ";
	}
	cout << endl;
}

// }}}

#endif

