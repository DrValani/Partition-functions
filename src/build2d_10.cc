#ifndef build2d_10_CC
#define build2d_10_CC

#include "latticeSize.cc"
#include "Storage.cc"
#include "Storage_y.cc"
#include "commonMethods.cc"
#include "searchMethods.cc"

int forward_traverse_tree(int tree_index, int current_column, Storage_y svy,
    vector<Storage> &config_tree, std::vector<int> &deltaCounter,
    const G_TABLE1 &ylist_bc3, const vector<Storage> & GrowthMatrix) {
	int configuration = config_tree[tree_index][svy[ySize - 1 - current_column]];
	if (current_column < ySize - 1) {
		return forward_traverse_tree(configuration, current_column + 1, svy,
		    config_tree, deltaCounter, ylist_bc3, GrowthMatrix);
	}

	if (configuration == -1) {
		svy.calc_subgraphs();
		deltaCounter.push_back(svy.getDelta(ylist_bc3, GrowthMatrix));
		config_tree[tree_index].push_back(svy[ySize - 1 - current_column],
		    deltaCounter.size() - 1);
	}
	return config_tree[tree_index][svy[ySize - 1 - current_column]];
}

void create_ylist_table2(int tree_index, const int current_column,
    G_TABLE2 &ylistTable, vector<int> yogi, vector<Storage> &config_tree) {
	if (current_column < ySize - 1) {
		for (int config = 0; config < ARRSIZE; config++) {
			if (config_tree[tree_index][config] != -1) {
				yogi[ySize - 1 - current_column] = config;
				create_ylist_table2(config_tree[tree_index][config], current_column + 1,
				    ylistTable, yogi, config_tree);
			}
		}
	} else {
		for (int config = 0; config < ARRSIZE; config++) {
			if (config_tree[tree_index][config] != -1) {
				yogi[ySize - 1 - current_column] = config;
				yogi[ySize] = config_tree[tree_index][config];
				ylistTable.push_back(yogi);
			}
		}
	}
}

/**
 *Difference is find_symmetrical_equivalent is used instead
 */
Storage lastColumn(G_TABLE1 &transformGroups, Storage_y svy,
    vector<Storage> &config_tree, std::vector<int> &deltaCounter,
    const G_TABLE1 &ylist_bc3, const vector<Storage> & GrowthMatrix) {
	Storage outer(true);
	for (int a = 0; a < ARRSIZE; a++) {
		int new_a = find_symmetrical_equivalent(a, transformGroups, outer);
		if (new_a == -1) {
			svy.push_back(0, a);
			Storage_y tempsvy = svy;
			int alternate = forward_traverse_tree(config_tree.size() - 1, 0, tempsvy,
			    config_tree, deltaCounter, ylist_bc3, GrowthMatrix);
			outer.push_back(a, alternate);
		} else outer.push_back(a, outer[new_a]);
	}
	return outer;
}

int buildBase(const int current_column, vector<Storage> &config_tree,
    G_TABLE1 transformGroups, vector<vector<Storage> > &vlist, Storage_y &svy,
    std::vector<int> &deltaCounter, vector<Storage> &vpol_vlist,
    const G_TABLE1 &ylist_bc3, const vector<Storage> & GrowthMatrix) {

	//if (current_column < ySize-1) {
	if (current_column < ySize - 1) {
		Storage souter(true);
		for (int config = 0; config < ARRSIZE; config++) {
			int sym_equivalent_config = find_symmetrical_equivalent(config,
			    transformGroups, souter);
			if (sym_equivalent_config == -1) {
				G_TABLE1 sub_transformGroups;
				Calculate_subgroup(sub_transformGroups, config, transformGroups);
				svy.push_back(ySize - 1 - current_column, config);
				int vector_index = buildBase(current_column + 1, config_tree,
				    sub_transformGroups, vlist, svy, deltaCounter, vpol_vlist,
				    ylist_bc3, GrowthMatrix);
				souter.push_back(config, vector_index);
			} else {
				souter.push_back(config, -1);
			}
		}
		vlist[0].push_back(souter);
		return vlist[0].size() - 1;
	} else {
		Storage s = lastColumn(transformGroups, svy, config_tree, deltaCounter,
		    ylist_bc3, GrowthMatrix);
		return search_add_list_exact(vpol_vlist, s);
	}
}

// {{{ build 2d lattice

/**
 * This is where the build starts.
 * Might have to send the vlist, grouplist objects back to main
 * I know it works for all xSizes >= 2.
 */
void build_2dLattice(vector<vector<Storage> > &vlist1, G_TABLE1 transformGroups,
    G_TABLE2 &ylistTable, vector<Polynomial> &vpol, vector<Storage> &vpol_vlist,
    const G_TABLE1 &ylist_bc3, const vector<Storage> &GrowthMatrix) {

	std::vector<int> deltaCounter;
	vector < Storage > config_tree;

	create_horizontal_configuration_tree(0, config_tree, transformGroups);

	Storage_y svy;
	buildBase(0, config_tree, transformGroups, vlist1, svy, deltaCounter,
	    vpol_vlist, ylist_bc3, GrowthMatrix);
	std::vector<int> yogi(ySize + 1, -1);

	create_ylist_table2(config_tree.size() - 1, 0, ylistTable, yogi, config_tree);

	make_vpol(vpol, deltaCounter);
	std::sort(ylistTable.begin(), ylistTable.end());
}
// }}}

#endif


