/* 
 * File:   Grid.h
 * Author: yogi
 *
 * Created on 17 August 2012, 07:11
 */
#ifndef GRID_H
#define	GRID_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include "latticeSize.cc"
#include "Globals.cc"
using namespace std;
/**
 * The binary class makes a Grid with n rows and m cols.  
 * It has a method add, which will go through             
 * every combination and return false once all possible   
 * combinations have been called.                         
 **/
class Grid{
public:  
  Grid():row_size(1),col_size(1){}
  Grid(const int &_row_size, const int &_col_size);
  bool equals(const Grid &g)const;
  void transform(const vector<int> &perm);
  void reflect();
  void shift();
  void rotate();
  bool combinationsRemaining() const;
  void print() const;
  // New methods
  int calcEnergy_step1()const;
  int calcEnergy_step2()const;
  int calcEnergy_step3()const;
  int calcEnergy_step3_1()const;
  int calcEnergy_step4()const;
  int calcEnergy_brute2d_1()const;
  int calcEnergy_brute2()const;
  int calcEnergy_brute3()const;
  int calcEnergy_brute4()const;
  int calcEnergy_brute5()const;
  int calcEnergy_brute6()const;
  int calcEnergy_brute7()const;
  int calcEnergy_brute8()const;
  int calcEnergy_brute9()const;
  void reset();
  void nextCombination();
  int getPoint(int x, int y)const{return (bv[x][y]);}
  int get_row_size()const{return row_size;}
  int get_col_size()const{return col_size;}
  void setGrid(const int &x,const int &y, const int &i){bv[x][y] = i;}
  //set point (x, y) to i
  void calculate_edge_subgraphs(vector<pair<int, int> > &vpr, bool inverse)const;
private:
  vector<vector<int> > bv; //This is the Grid
  int row_size;
  int col_size;
  bool moreCombos;  //Are there any more combinations?
  void changeAll(const int change, const int to)const;
  void check_neighbour(const int &i,const int &j, const int &ngh_i,
		       const int &ngh_j, int &subG_counter, int &subG_no,
		       vector<pair<int, int> > &vpr, bool inverse)const;
}; //End of Class Grid 
#endif	/* GRID_H */

