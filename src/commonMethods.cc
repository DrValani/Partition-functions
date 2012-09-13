#ifndef COMMONMETHODS_CC
#define COMMONMETHODS_CC
#include <vector>
#include "Grid.cc"
#include "Storage.cc"
#include "Polynomial.cc"
#include "refPolynomial.cc"
using namespace std;

void final_check(Polynomial &p, const int &z) {
    mpz_class total = 0;
    int MAXPOLYSIZE = (3 * xSize * ySize * zSize) + 1;
    for (int i = MAXPOLYSIZE - 1; i >= 0; i--) {
        total = total + p[i];
    }

    mpz_class total2;
    mpz_ui_pow_ui(total2.get_mpz_t(), Q, xSize * ySize * z);
    total = total - total2;
    if (total != 0) {
        cout << "total is: " << total << endl;
        cout << "The total does not match the total number of combos\n QBAR = " << ARRSIZE << "\n";
        throw new runtime_error("The total does not match the total number of combos\n");
    }
}

void make_vpol(vector<Polynomial> &vpol, vector<int> &ylist) {
    vpol.clear();
	int MAXPOLYSIZE =(3 * xSize * ySize * zSize) + 1;
	Polynomial p(MAXPOLYSIZE);
    for (int i = 0; i < ylist.size(); i++) {
        p.push_back(1, ylist[i]);
        vpol.push_back(p);
        p.clear();
    }
}

void check_ylist_table(G_TABLE2 &ylistTable) {
    cout << "now checking ylist\n";
    for (int i = 0; i < ylistTable.size(); i++)
        if (ylistTable[i][ySize] != i) {
            cout << "ERROR: check_ylist_table\n";
            cout << "cannot use ylist[i][ySize] as polynomial multiplicity\n";
            throw new runtime_error("ERROR: check_ylist_table");
        }
    cout << "done\n";
}

void Create_energy_list(vector<Storage> &elist) {
    Grid y(xSize, 2);
    int delta;
    for (int i = 0; i < ARRSIZE; i++) {
        Storage s;
        for (int j = 0; j < ARRSIZE; j++) {
            delta = y.calcEnergy_step3();
            s.push_back(j, delta);
            y.nextCombination();
        }
        elist.push_back(s);
    }
}

/**
 * This method creates the first initial line, 
 *
 *   0
 *   0
 *   1
 *   
 *  all delta values are calculated and stored in vlist
 *
 **/

void Create_energy_list(Storage &elist) {
    Grid y(xSize, 1);
    int delta;
    for (int j = 0; j < ARRSIZE; j++) {
        delta = y.calcEnergy_step3_1();
        elist.push_back(j, delta);
        y.nextCombination();
    }
}

void Create_energy_matrix(vector<Storage>   &GrowthMatrix) {
    Create_energy_list(GrowthMatrix);
    Storage s;
    Create_energy_list(s);
    GrowthMatrix.push_back(s);
}

void Calculate_subgroup(G_TABLE1 &sub_transformGroups, const int &index,
        const G_TABLE1 &transformGroups) {
    for (int pg = 0; pg < transformGroups.size(); pg++) {
        int gn = TRANSFORMGROUP[transformGroups[pg]][index];
        if (gn == index && transformGroups[pg] != 0) {
            sub_transformGroups.push_back(TRANSFORMGROUP[transformGroups[pg]].getId());
        }
    }
}

/**
 *This method find horizontal symmetrical equivalent configurations.
 *@returns -1 if there is no symmertrical equivalent or that configuration has not been found yet
 * Eg for c = 0, and c' = 31, it will return -1 for c and 0 for c'.
 */
int find_symmetrical_equivalent(const int &index, const G_TABLE1 &transformGroups,
        const Storage &outer) {
    for (int transformGroup_index = 0; transformGroup_index < transformGroups.size(); transformGroup_index++) {
        int symmetrical_equivalent = TRANSFORMGROUP[transformGroups[transformGroup_index]][index];
        if (outer[symmetrical_equivalent] != -1) {
            return symmetrical_equivalent;
        }
    }
    return -1;
}

/**
 *This method is the inverse of find_symmetrical_equivalent configurations.
 *@returns The symmetrical equivalent form of the require configuration
 *@param also transformGroup is converts c' to c is set.
 * Eg for c = 0, and c' = 31, it will return -1 for c and 0 for c'.
 */
int inverse_find_symmetrical_equivalent(const int &config, const G_TABLE1 &transformGroups,
        const Storage &outer, int &transformGroupOperation) {
    for (int transformGroup_index = 0; transformGroup_index < transformGroups.size(); transformGroup_index++) {
        int symmetrical_equivalent = TRANSFORMGROUP[transformGroups[transformGroup_index]][config];
        if (outer[symmetrical_equivalent] != -1) {
            // Found it, config is equivalent to symmetrical_eqivalent
            transformGroupOperation = TRANSFORMGROUP[transformGroup_index].getId();
            return symmetrical_equivalent;
        }
    }
    return -1;
}

/**
 *This method will replace the hash table. Will create a tree, for configuration
 * c = (s1,...sy). Where the yth node will be a reference to its corresponding polynomial.
 *
 */
int create_horizontal_configuration_tree(const int current_column, vector<Storage>   &config_tree,
        G_TABLE1 transformGroups) {

    if (current_column < ySize - 1) {
        Storage souter(true);
        for (int config = 0; config < ARRSIZE; config++) {
            int sym_equivalent_config = find_symmetrical_equivalent(config, transformGroups, souter);
            if (sym_equivalent_config == -1) {
                G_TABLE1 sub_transformGroups;
                Calculate_subgroup(sub_transformGroups, config, transformGroups);
                int vector_index = create_horizontal_configuration_tree(
                        current_column + 1, config_tree, sub_transformGroups);
                souter.push_back(config, vector_index);
            } else {
                souter.push_back(config, -1);
            }
        }
        config_tree.push_back(souter);
        return config_tree.size() - 1;
    } else {
        Storage s(true);
        config_tree.push_back(s);
        return config_tree.size() - 1;
    }
}

int create_horizontal_configuration_tree2(const int current_column, vector<Storage>   &config_tree,
        G_TABLE1 transformGroups) {

    if (current_column > 0) {
        Storage souter(true);
        for (int config = 0; config < ARRSIZE; config++) {
            int sym_equivalent_config = find_symmetrical_equivalent(config, transformGroups, souter);
            if (sym_equivalent_config == -1) {
                G_TABLE1 sub_transformGroups;
                Calculate_subgroup(sub_transformGroups, config, transformGroups);
                int vector_index = create_horizontal_configuration_tree(
                        current_column - 1, config_tree, sub_transformGroups);
                souter.push_back(config, vector_index);
            } else {
                souter.push_back(config, -1);
            }
        }
        config_tree.push_back(souter);
        return config_tree.size() - 1;
    } else {
        Storage s(true);
        config_tree.push_back(s);
        return config_tree.size() - 1;
    }
}

/**
 *This method sums up the vlist
 *
 */
int make_refpol(
        const int current_column,
        G_TABLE1 transformGroups,
        vector<refPolynomial> &vrefpol,
        vector<Storage>   &vlist,
        vector<Storage>   &vpol_vlist,
        int vlist_index
        ) {


    vector<refPolynomial> vrefpolnew;
    refPolynomial refpol;
    if (current_column < ySize - 1) {
        Storage souter(true);
        for (int config = 0; config < ARRSIZE; config++) {
            if (vlist[vlist_index][config] != -1) {
                G_TABLE1 sub_transformGroups;
                Calculate_subgroup(sub_transformGroups, config, transformGroups);

                int vector_index = make_refpol(
                        current_column + 1,
                        sub_transformGroups,
                        vrefpolnew,
                        vlist,
                        vpol_vlist,
                        vlist[vlist_index][config]
                        );
                souter.push_back(config, vector_index);
                refpol.push_back(vector_index, 1);

            } else {
                int sym_equivalent_config = find_symmetrical_equivalent(config, transformGroups, vlist[vlist_index]);
                refpol.push_back(souter[sym_equivalent_config], 1);
            }
        }
    } else {
        for (int config = 0; config < ARRSIZE; config++)
            refpol.push_back(vpol_vlist[vlist_index][config], 1);
        vrefpol.push_back(refpol);
        refpol.clear();
    }
    if (current_column < ySize - 1) {
        refPolynomial refpolnew;
        refpolnew.sum(refpol, vrefpolnew);
        vrefpol.push_back(refpolnew);
    }
    return vrefpol.size() - 1;
}

void clearLists(vector<vector<Storage> >  &vvlist) {
    for (vector<vector<Storage> > ::iterator i = vvlist.begin(); i != vvlist.end(); ++i)
        i->clear();
}

void clearLists(vector<vector<Polynomial> > &vvlist) {
    for (vector<vector<Polynomial> >::iterator i = vvlist.begin(); i != vvlist.end(); ++i)
        i->clear();
}

int find_ylist_reference(const G_TABLE2 &ylist, const G_TABLE1 &ylist_bc3) {
    for (int i = 0; i < ylist.size(); i++) {
        bool equal = true;
        for (int j = 0; j < ySize; j++)
            if (ylist[i][j] != ylist_bc3[j]) {
                equal = false;
                break;
            }
        if (equal) {
            return ylist[i][ySize];
        }
    }
    throw new runtime_error("Error:z_pol_ref\n");
}

int calculate_outer_hamiltonian(const vector<Storage>   &GrowthMatrix, const G_TABLE1 &current_config) {
    int Hamiltonian = 0;
    for (int i = 0; i < ySize; i++) {
        Hamiltonian += GrowthMatrix[ARRSIZE][current_config[i]]; // this is the line
    }
    for (int i = 1; i < ySize; i++) {
        Hamiltonian += GrowthMatrix[current_config[i]][current_config[i - 1]];
    }
    if (BConds1)
        Hamiltonian += GrowthMatrix[current_config[0]][current_config[ySize - 1]];
    return Hamiltonian;
}

void addEnergy4pol(const Storage &s, const Storage &e,
        refPolynomial &refpol, const int &maxpower) {
    int expo;
    refpol.clear();
    for (int i = 0; i < ARRSIZE; i++) {
        expo = s[i] * maxpower + e[i];
        refpol.push_back(expo, 1);
    }
}
int SMALLEST_TRANSFERMATRIX = 0;

void calculate_inner_vpol(const int &current_column,
        vector<vector<Polynomial> > &vpol,
        vector<vector<refPolynomial> >  &vrefpol,
        Polynomial &p
        ) {
    if (vpol[current_column + 2].size() < vrefpol[current_column].size()) {
        p.clear();
        while (vpol[current_column + 2].size() < vrefpol[current_column].size()) {
            vpol[current_column + 2].push_back(p);
        }
    }
    vector<vector<refPolynomial > > vrefpolnew(ySize + 1);
    refPolynomial refpol;
    int mxpow = 1;
    int cc = 0;
    int starter = 0;
    if (SMALLEST_TRANSFERMATRIX < current_column)
        cc = SMALLEST_TRANSFERMATRIX;
lookagain2:
    while (cc < current_column && vrefpol[cc].size() == 0)
        cc++;
    int cc1 = cc + 1;
    while (cc1 < current_column && vrefpol[cc1].size() == 0) {
        cc = cc1;
        goto lookagain2;
    }
    bool SWAP = true;
    if (cc < current_column && vrefpol[cc].size() > 0) {
        vrefpolnew[cc].swap(vrefpol[cc]);
        starter = cc;
        while (cc < current_column) {
            for (int i = 0; i < vrefpol[cc + 1].size(); i++) {
                refpol.clear();
                refpol.sum(vrefpol[cc + 1][i], vrefpolnew[cc], MAXPOWER*mxpow, MAXPOWER);
                vrefpolnew[cc + 1].push_back(refpol);
            }
            if (SWAP)
                vrefpolnew[cc].swap(vrefpol[cc]);
            SWAP = false;
            mxpow++;
            cc++;
        }
    } else {
        vrefpolnew[current_column].swap(vrefpol[current_column]);
        starter = current_column;
    }

    p.clear();
    for (int i = 0; i < vrefpolnew[current_column].size(); i++) {
        vpol[current_column + 2][i].clear();
        vrefpolnew[current_column][i].sum(vpol[current_column + 2][i], vpol[starter + 1], MAXPOWER * mxpow, p);
    }
}
#endif
