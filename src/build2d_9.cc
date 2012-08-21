#ifndef build2d_9_CC
#define build2d_9_CC

#include "latticeSize.cc"
#include "commonMethods.cc"
#include "searchMethods.cc"
#include "myTimer.cc"
#include "Storage_y.cc"

int forward_traverse_tree(int tree_index, int current_column, Storage_y svy,
    vector<Storage> &config_tree, vector<int> &deltaCounter) {
	int configuration = config_tree[tree_index][svy[current_column]];
	if (current_column < ySize - 1) {
		return forward_traverse_tree(configuration, current_column + 1, svy,
		    config_tree, deltaCounter);
	} else if (configuration == -1) {
		svy.calc_subgraphs();
		deltaCounter.push_back(svy.getDelta());
		config_tree[tree_index].push_back(svy[current_column],
		    deltaCounter.size() - 1);
	}
	return config_tree[tree_index][svy[current_column]];
}

void create_ylist_table(int tree_index, const int current_column,
    G_TABLE2 &ylistTable, vector<int> yogi, vector<Storage> &config_tree) {
	if (current_column < ySize - 1) {
		for (int config = 0; config < ARRSIZE; config++) {
			if (config_tree[tree_index][config] != -1) {
				yogi[current_column] = config;
				create_ylist_table(config_tree[tree_index][config], current_column + 1,
				    ylistTable, yogi, config_tree);
			}
		}
	} else {
		for (int config = 0; config < ARRSIZE; config++) {
			if (config_tree[tree_index][config] != -1) {
				yogi[current_column] = config;
				yogi[ySize] = config_tree[tree_index][config];
				ylistTable.push_back(yogi);
				for (int i = 0; i <= ySize; i++)
					cout << yogi[i] << " ";
				cout << endl;
			}
		}
	}
}

/**
 *Difference is find_symmetrical_equivalent is used instead
 */
Storage lastColumn(G_TABLE1 &transformGroups, Storage_y svy,
    vector<Storage> &config_tree, vector<int> &deltaCounter) {
	Storage outer(true);
	for (int a = 0; a < ARRSIZE; a++) {
		int new_a = find_symmetrical_equivalent(a, transformGroups, outer);
		if (new_a == -1) {
			svy.push_back(0, a);
			Storage_y tempsvy = svy;
			tempsvy.calc_subgraphs();
			tempsvy.setSmallest();

			int alternate = forward_traverse_tree(config_tree.size() - 1, 0, tempsvy,
			    config_tree, deltaCounter);
			outer.push_back(a, alternate);
		} else outer.push_back(a, outer[new_a]);
	}
	return outer;
}

int buildBase(const int current_column, vector<Storage> &config_tree,
    G_TABLE1 transformGroups, vector<vector<Storage> > &vlist, Storage_y &svy,
    vector<int> &deltaCounter, vector<Storage> &vpol_vlist) {

	if (current_column > 0) {
		Storage souter(true);
		for (int config = 0; config < ARRSIZE; config++) {

			int sym_equivalent_config = find_symmetrical_equivalent(config,
			    transformGroups, souter);
			if (sym_equivalent_config == -1) {
				G_TABLE1 sub_transformGroups;
				svy.push_back(current_column, config);
				Calculate_subgroup(sub_transformGroups, config, transformGroups);
				int vector_index = buildBase(current_column - 1, config_tree,
				    sub_transformGroups, vlist, svy, deltaCounter, vpol_vlist);
				souter.push_back(config, vector_index);
			} else {
				souter.push_back(config, -1);
			}
			if (current_column == 0) calculateETA(startTime, config + 1, ARRSIZE);
		}
		vlist[0].push_back(souter);
		return vlist[0].size() - 1;
	} else {
		Storage s = lastColumn(transformGroups, svy, config_tree, deltaCounter);
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
    G_TABLE2 &ylistTable, vector<Polynomial> &vpol,
    vector<Storage> &vpol_vlist) {

	vector<int> deltaCounter;
	vector < Storage > config_tree;
	time (&startTime);
	cout << " create horizontal tree " << endl;
	create_horizontal_configuration_tree(0, config_tree, transformGroups);
	Storage_y svy;
	time(&startTime);
	cout << " create base layer" << endl;
	buildBase(ySize - 1, config_tree, transformGroups, vlist1, svy, deltaCounter,
	    vpol_vlist);
	vector<int> yogi(ySize + 1);
	time(&startTime);
	cout << " create ylist Table" << endl;
	create_ylist_table(config_tree.size() - 1, 0, ylistTable, yogi, config_tree);
	cout << "yist.size() " << deltaCounter.size() << endl;

	make_vpol(vpol, deltaCounter);
	cout << "now sorting ylist" << endl;
	sort(ylistTable.begin(), ylistTable.end());
	check_ylist_table(ylistTable);
}
// }}}

#endif

