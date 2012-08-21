#ifndef buildGroups_CC
#define buildGroups_CC

#include "Globals.cc"
#include "Storage.cc"
#include "Grid.cc"
#include "searchMethods.cc"
// {{{ makeGroups
void transformlist(Storage o, Storage &r, const int &i) {
	for (int j = 0; j < ARRSIZE; j++)
		r.push_back(j, o[transformARR[i][j]]);
}
void reflectlist(Storage o, Storage &r) {
	for (int i = 0; i < ARRSIZE; i++)
		r.push_back(i, o[reflectxARR[i]]);
}
void shiftlist(Storage o, Storage &r) {
	for (int i = 0; i < ARRSIZE; i++)
		r.push_back(i, o[shiftxARR[i]]);
}
void shiftAll(Storage &identity, Storage shift) {

	for (int i = 0; i < xSize; i++) {
		Storage g4;
		g4.createGroup(identity, shift);
		g4.searchGroup(TRANSFORMGROUP);
		Storage transform;
		for (int j = 0; j < transformARR.size(); j++) {
			transformlist(shift, transform, j);
			Storage g5;
			g5.createGroup(identity, transform);  // change 0 to 1, 1 to 0
			g5.searchGroup(TRANSFORMGROUP);
		}
		Storage reflect;
		reflectlist(shift, reflect);
		Storage g6;
		g6.createGroup(identity, reflect);  // reflect in x axis
		g6.searchGroup(TRANSFORMGROUP);
		Storage rT;
		for (int j = 0; j < transformARR.size(); j++) {
			transformlist(reflect, rT, j);
			Storage g7;
			g7.createGroup(identity, rT);  // reflect in x axis and transform
			g7.searchGroup(TRANSFORMGROUP);
		}
		shiftlist(shift, shift);
	}
}
/**
 *Creates an n x n grouptable of all the groups
 **/
void makeGroupTable() {
	int numberOfGroups = TRANSFORMGROUP.size();
	if (true) {
		G_TABLE1 rows(numberOfGroups, -1);
		for (int i = 0; i < numberOfGroups; i++) {
			GROUPTABLE.push_back(rows);
			TRANSFORMGROUP[i].setId(i);
		}
	}
	Storage identity;
	for (int i = 0; i < numberOfGroups; i++) {
		for (int j = 0; j < numberOfGroups; j++) {
			identity = TRANSFORMGROUP[i].applyGroup(TRANSFORMGROUP[j]);
			int found = search_add_list_exact(TRANSFORMGROUP, identity);
			GROUPTABLE[i][j] = found;
			cout << found << " ";
		}
		cout << endl;
	}
	if (TRANSFORMGROUP.size() != numberOfGroups) {
		cout << "Not all groups had been found. Review method for find groups\n";
		exit(0);
	}
}
;

// }}}

// {{{ makeGroup reflect shift transpose

void makeGroup_reflectx() {
	// makeGroup_transpose has to be run before this function
	Grid g(xSize, 1);
	int gindex = 0;
	while (g.combinationsRemaining()) {
		g.reflect();
		bool found = false;
		int i = 0;
		while (!found) {
			bool stop = false;
			for (int j = 0; j < xSize && !stop; j++) {
				if (g.getPoint(j, 0) != transposeARR[i][j]) stop = true;
			}
			if (!stop) {
				reflectxARR[gindex] = i;
				found = true;
			} else i++;
		}
		g.reflect();
		g.nextCombination();
		gindex++;
	}
}
void makeGroup_shiftx() {
	// makeGroup_shiftx has to be run before this function
	if (!BConds2)  // If the boundary condition is false then a shift will have no
	//effect. So make shiftxARR the identity
	for (int i = 0; i < ARRSIZE; i++)
		shiftxARR[i] = i;
	else {
		Grid g(xSize, 1);
		Grid g_temp(xSize, 1);
		int gindex = 0;
		while (g.combinationsRemaining()) {
			g_temp = g;
			g_temp.shift();
			bool found = false;
			int i = 0;
			while (!found) {
				bool stop = false;
				for (int j = 0; j < xSize && !stop; j++) {
					if (g_temp.getPoint(j, 0) != transposeARR[i][j]) stop = true;
				}
				if (!stop) {
					shiftxARR[gindex] = i;
					found = true;
				} else i++;
			}
			g.nextCombination();
			gindex++;
		}
	}
}
void makeGroup_transpose() {
	Grid g(xSize, 1);
	for (int i = 0; i < ARRSIZE; i++) {
		for (int j = 0; j < xSize; j++)
			transposeARR[i][j] = g.getPoint(j, 0);
		g.nextCombination();
	}
}
void makeGroup_transform() {
	vector<int> tes;
	for (int i = 0; i < Q; i++)
		tes.push_back(i);
	while (next_permutation(tes.begin(), tes.end())) {
		Grid g(xSize, 1);
		Grid g_temp(xSize, 1);
		vector<int> transGroup;
		while (g.combinationsRemaining()) {
			g_temp = g;
			g_temp.transform(tes);
			bool found = false;
			int i = 0;
			while (!found) {
				bool stop = false;
				for (int j = 0; j < xSize && !stop; j++) {
					if (g_temp.getPoint(j, 0) != transposeARR[i][j]) stop = true;
				}
				if (!stop) {
					transGroup.push_back(i);
					found = true;
				} else i++;
			}
			g.nextCombination();
		}
		transformARR.push_back(transGroup);
		transGroup.clear();
		cout << endl;
		for (int i = 0; i < Q; i++)
			cout << tes[i] << " ";
		cout << endl;
	}
}

// }}}
void makeGroups() {
	makeGroup_transpose();
	makeGroup_reflectx();
	makeGroup_shiftx();
	makeGroup_transform();

	Storage identity;
	for (int i = 0; i < ARRSIZE; i++)
		identity.push_back(i, i);
	Storage g;
	g.createGroup(identity, identity);  // identity e
	g.searchGroup(TRANSFORMGROUP);

	Storage transform;
	for (int j = 0; j < transformARR.size(); j++) {
		transformlist(identity, transform, j);
		Storage g1;
		g1.createGroup(identity, transform);  // change 0 to 1, 1 to 0
		g1.searchGroup(TRANSFORMGROUP);
	}

	Storage reflect;
	reflectlist(identity, reflect);
	Storage g2;
	g2.createGroup(identity, reflect);  // reflect in x axis
	g2.searchGroup(TRANSFORMGROUP);

	Storage rT;
	for (int j = 0; j < transformARR.size(); j++) {
		transformlist(reflect, rT, j);
		Storage g3;
		g3.createGroup(identity, rT);  // reflect in x axis and transform
		g3.searchGroup(TRANSFORMGROUP);
	}
	if (BConds2) {
		Storage shift = identity;
		for (int i = 0; i < xSize; i++) {
			shiftAll(identity, shift);
			for (int j = 0; j < transformARR.size(); j++) {
				transformlist(shift, transform, j);
				Storage g1;
				g1.createGroup(identity, transform);  // change 0 to 1, 1 to 0
				g1.searchGroup(TRANSFORMGROUP);
				shiftAll(identity, transform);
			}
			reflectlist(shift, reflect);
			Storage g2;
			g2.createGroup(identity, reflect);  // reflect in x axis
			g2.searchGroup(TRANSFORMGROUP);
			shiftAll(identity, reflect);

			for (int j = 0; j < transformARR.size(); j++) {
				transformlist(reflect, rT, j);
				Storage g3;
				g3.createGroup(identity, rT);  // reflect in x axis and transform
				g3.searchGroup(TRANSFORMGROUP);
				shiftAll(identity, rT);
			}
			shiftlist(shift, shift);
		}
	}
	makeGroupTable();

}

#endif

