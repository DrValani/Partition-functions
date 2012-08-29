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
  bool combinationsRemaining() const;
  void print() const;
  // New methods
  int calcEnergy_step3()const;
  int calcEnergy_step3_1()const;
  int calcEnergy_brute2d_1()const;
  void reset();
  void nextCombination();
  int getPoint(int x, int y)const{return (bv[x][y]);}
  void setGrid(const int &x,const int &y, const int &i){bv[x][y] = i;}
  int get_row_size()const{return row_size;}
  int get_col_size()const{return col_size;}
  //set point (x, y) to i
private:
  vector<vector<int> > bv; //This is the Grid
  int row_size;
  int col_size;
  bool moreCombos;  //Are there any more combinations?
}; //End of Class Grid 
#endif	/* GRID_H */

