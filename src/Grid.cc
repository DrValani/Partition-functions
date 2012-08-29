#ifndef GRID_CC
#define GRID_CC

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#include "Grid.h"

// {{{ calcEnergy

int Grid::calcEnergy_step3() const {
	int delta = 0;
	for (int indexR = 0; indexR < row_size; indexR++) {
		if (bv[indexR][0] == bv[indexR][1]) {
			delta++;
		}
	}
	return delta;
}
int Grid::calcEnergy_step3_1() const {
	int delta = 0;
	for (int indexR = 0; indexR < row_size - 1; indexR++) {
		if (bv[indexR][0] == bv[indexR + 1][0]) {
			delta++;
		}
	}
	if (BConds2)
		if (bv[row_size - 1][0] == bv[0][0]) {
			delta++;
		} // (1)
	return delta;
}

int Grid::calcEnergy_brute2d_1() const {
	int delta = 0;
	for (int indexC = 0; indexC < col_size; indexC++) {
		for (int indexR = 0; indexR < row_size - 1; indexR++) {
			if (bv[indexR][indexC] == bv[indexR + 1][indexC]) {
				delta++;
			}
		}
	}
	for (int indexR = 0; indexR < row_size; indexR++) {
		for (int indexC = 0; indexC < col_size - 1; indexC++) {
			if (bv[indexR][indexC] == bv[indexR][indexC + 1]) {
				delta++;
			}
		}
	}
	if (BConds2)
		for (int indexC = 0; indexC < col_size; indexC++) {
			if (bv[row_size - 1][indexC] == bv[0][indexC]) {
				delta++;
			} // (1)
		}
	if (BConds1)
		for (int indexR = 0; indexR < row_size; indexR++) {
			if (bv[indexR][0] == bv[indexR][col_size - 1]) {
				delta++;
			}
		}
	return delta;
}
// }}}

// {{{ Methods

bool Grid::equals(const Grid &g) const {
	for (int r = 0; r < row_size; r++) {
		for (int c = 0; c < col_size; c++) {
			if (this->bv[r][c] != g.bv[r][c]) {
				return false;
			}
		}
	}
	return true;
}
void Grid::transform(const vector<int> &perm) {
	for (int r = 0; r < row_size; r++)
		for (int c = 0; c < col_size; c++)
			this->bv[r][c] = perm[this->bv[r][c]];
}
void Grid::reflect() {
	const Grid temp = *this;
	for (int r = 0; r < row_size; r++) {
		for (int c = 0; c < col_size; c++) {
			this->bv[r][c] = temp.bv[row_size - 1 - r][c];
		}
	}
}
void Grid::shift() {
	const Grid temp = *this;
	for (int c = 0; c < col_size; c++) {
		for (int r = 0; r < row_size - 1; r++) {
			this->bv[r + 1][c] = temp.bv[r][c];
		}
		this->bv[0][c] = temp.bv[row_size - 1][c];
	}
}
void Grid::reset() {
	for (int indexR = 0; indexR < row_size; indexR++) {
		for (int indexC = 0; indexC < col_size; indexC++) {
			bv[indexR][indexC] = 0;
		}
	}
	moreCombos = true;
}

/**********************************************\
 *   Change the Grid to the next combination. *
 *   Increses from top to bottom,             *
 *   then left to right.                      *
 *   So the last combination is on the        *
 *   bottom right corner of the Grid.         *
 *                                            *
 *   returns false if all possible            *
 *   combinations have been called            *
 \**********************************************/

void Grid::nextCombination() {
	int indexR = 0;
	int indexC = 0;
	bool stop = false;

	while (!stop) {
		if (bv[indexR][indexC] < Q - 1) {
			bv[indexR][indexC]++;
			stop = true;
		} else {
			bv[indexR][indexC] = 0;
			if (indexR < row_size - 1) {
				indexR++;
			} else if (indexC < col_size - 1) {
				indexC++;
				indexR = 0;
			} else { //No more combinations
				stop = true;
				moreCombos = false;
			}
		}
	}
}

/* Make a n x m Grid. */
Grid::Grid(const int &_row_size, const int &_col_size) :
		row_size(_row_size), col_size(_col_size), moreCombos(true) {

	for (int indexR = 0; indexR < _row_size; indexR++) {
		vector<int> temp;
		for (int indexC = 0; indexC < _col_size; indexC++) {
			temp.push_back(0);
		}
		bv.push_back(temp);
	}
}

/*Is there any more combinations left, return true if there are. */
bool Grid::combinationsRemaining() const {
	return moreCombos;
}

/* Print out the Grid. */
void Grid::print() const {
	for (int indexR = 0; indexR < row_size; indexR++) {
		for (int indexC = 0; indexC < col_size; indexC++) {
			cout << bv[indexR][indexC];
		}
		cout << endl;
	}
	cout << endl;
}

// {{{ Test class method

/*Test out the class */
/*int main(){
 cout << "start" << endl;
 float x = 0; float power = 0;
 Grid Fourth;
 while(Fourth.combinationsRemaining()){
 Grid Third;
 while(Third.combinationsRemaining()){
 Grid Second;
 while(Second.combinationsRemaining()){
 Grid First;
 while(First.combinationsRemaining()) {
 power = 0;
 power = First.calcEnergyClosed()
 + First.calcEnergyClosed(Second)
 + Second.calcEnergyClosed();

 cout << "Value: " << First.calcEnergyClosed()
 + First.calcEnergyClosed(Second)
 + Second.calcEnergyClosed() <<endl;

 First.nextCombination();
 x++;
 if(x == 262143){
 First.print();   Second.print();
 cout <<" AHHH" << endl;   }
 }
 Second.nextCombination();
 }
 Third.nextCombination();
 }
 Fourth.print();
 Fourth.nextCombination();
 }
 cout << "highest power: " << power << endl;
 cout << "total combos: " << x << endl;

 system("PAUSE");
 return 0;

 }*/

// }}}
// }}}
ostream & operator<<(ostream & out, const Grid &g) {
	int row_size = g.get_row_size();
	int col_size = g.get_col_size();
	for (int r = 0; r < row_size; r++) {
		for (int c = 0; c < col_size; c++) {
			out << g.getPoint(r, c);
		}
		out << endl;
	}
	return out;
}

ostream & operator<<(ostream & out, const vector<Grid> &g) {
	int row_size = g[0].get_row_size();
	int col_size = g[0].get_col_size();
	for (int r = 0; r < row_size; r++) {
		for (vector<Grid>::size_type s = 0; s < g.size(); s++) {
			for (int c = 0; c < col_size; c++) {
				out << g[s].getPoint(r, c);
			}
			out << " ";
		}
		out << endl;
	}
	return out;
}

// }}}
#endif
