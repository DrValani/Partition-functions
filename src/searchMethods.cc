#ifndef SEARCHMETHODS_CC
#define SEARCHMETHODS_CC
#include "Globals.cc"

template <class T, class S>
int search_list(T &elist, S &s) {
    int i = elist.size() - 1;
    while (i >= 0) {
        //  while(i < elist.size()){
        if (s == elist[i]) {
            return i;
        } else
            i--;
    }
    return -1;
}/**
 * Search the vector energy Lists if found return its place
 * else return -1;
 */

template <class T, class S>
int search_list_exact(T &elist, S &s) {
    int i = elist.size() - 1;
    while (i >= 0) {
        //  while(i < elist.size()){
        if (s.equals(elist[i])) {
            return i;
        } else
            i--;
    }
    return -1;
}

/**
 * Binary search. vs2 must be sorted.
 * returns an interger representing where it is found
 * if not found inserts into correct place;
 */

/*
   Search the vector energy Lists for a list s
   else add to list
   return where the storage s is stored
 */
template <class T, class S>
int search_add_list(T &elist, S &s) {
    int i = elist.size() - 1;
    //  while(i < elist.size()){
    while (i >= 0) {
        if (s == elist[i])
            return i;
        else
            i--;
    }
    elist.push_back(s);
    return elist.size() - 1;
}

template <class T, class S>
int search_add_list_exact(T &elist, S &s) {
    int i = elist.size() - 1;
    //  while(i < elist.size()){
    while (i >= 0) {
        if (s.equals(elist[i])) {
            return i;
        } else
            i--;
    }
    elist.push_back(s);
    return elist.size() - 1;
}

/**
 *  Search the vector of energy Lists for s.
 *  it will find the jth one
 *  if found return its place
 *  else return -1;
 */
template <class T, class S>
int search_add_list(T &elist, S &s, G_TABLE1 &hst2) {
    for (G_TABLE1::iterator i = hst2.begin(); i < hst2.end(); ++i) {
        if (*i != -1)
            if (s == elist[*i])
                return *i;
    }
    elist.push_back(s);
    return elist.size() - 1;
}
#endif
