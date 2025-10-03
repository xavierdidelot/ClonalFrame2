//
// Created by Didelot, Xavier on 29/09/2025.
//

#include "Tree.h"
#include "Alignment.h"

#include <iostream>

Tree::Tree(int n) {
    this->n = n;
    this->seq=nullptr;
    nodes = vector<Node *>();
    for (int i = 0; i < n; i++) {
        Node *node = new Node();
        node->set_name("seq" + std::to_string(i + 1));
        node->set_dist_father(1.0);
        nodes.push_back(node);
    }
    from_root = vector<Node *>();
    from_root.push_back(nodes[0]);
    from_root.push_back(nodes[1]);
    from_root.push_back(nodes[n - 1]);

    for (int i = n - 2; i >= 2; i--) {
        Node *node = new Node();
        nodes.push_back(node);
        node->set_left(nodes[i]);
        nodes[i]->set_father(node);
        node->set_dist_father(1.0);
        node->set_right(from_root[2]);
        from_root[2]->set_father(node);
        from_root[2] = node;
    }
}

string Tree::to_string() {
    return ("(" + from_root[0]->to_string() + ":" + std::to_string(from_root[0]->get_dist_father()) + "," +
            from_root[1]->to_string() + ":" + std::to_string(from_root[1]->get_dist_father()) + "," +
            from_root[2]->to_string() + ":" + std::to_string(from_root[2]->get_dist_father()) + ");");
}

Tree::~Tree() {
}

int Tree::get_n() {
    return n;
}

Node *Tree::get_node(int i) {
    return nodes[i];
}

void Tree::write(string filename) {
    ofstream file;
    file.open(filename);
    if (file.is_open() == false) {
        cout << "Error opening file " << filename << endl;
        exit(1);
    }
    file << to_string() << endl;
    file.close();
}

bool Tree::check_structure() {
    if (this->n * 2 - 3 - nodes.size() != 0) return false;
    if (from_root.size() != 3) return false;
    for (int i = 0; i < 3; i++) if (from_root[i]->get_father() != nullptr) return false;
    for (unsigned int i = 0; i < nodes.size(); i++) {
        if (nodes[i]->get_left() != nullptr && nodes[i]->get_left()->get_father() != nodes[i]) return false;
        if (nodes[i]->get_right() != nullptr && nodes[i]->get_right()->get_father() != nodes[i]) return false;
    }
    return true;
}

void Tree::add_mutations(int l, double theta) {
    seq = new string(l, 'A');
    uniform_int_distribution unif4(0, 3);
    for (int i = 0; i < l; i++) {
        *(seq->data() + i) = *(lets + unif4(rng));
    }
    for (int i = 0; i < 3; i++) from_root[i]->add_mutations(seq, theta);
}

void Tree::randomise(double avg_bra_len) {
    vector<Node *> to_merge = vector<Node *>();
    exponential_distribution rexp(avg_bra_len);
    for (int i = 0; i < n; i++) to_merge.push_back(nodes[i]);
    int m=n;
    while (to_merge.size() > 3) {
        uniform_int_distribution unif(0, (int) (to_merge.size() - 1));
        int a = unif(rng);
        int b = a;
        while (a == b) b = unif(rng);
        nodes[m]->set_left(to_merge[a]);
        nodes[m]->set_right(to_merge[b]);
        to_merge[a]->set_father(nodes[m]);
        to_merge[b]->set_father(nodes[m]);
        nodes[m]->set_dist_father(rexp(rng));
        to_merge[a] = nodes[m];
        to_merge[b] = to_merge.back();
        to_merge.pop_back();
        m++;
    }
    from_root=vector<Node *>();
    for (int i = 0; i < 3; i++) {
        to_merge[i]->set_father(nullptr);
        from_root.push_back(to_merge[i]);
    }
}

void Tree::nj(Alignment *ali) {
    vector<int> to_merge = vector<int>();
    for (int i = 0; i < n; i++) to_merge.push_back(i);
    vector<vector<double> > dist(2 * n, vector<double>(2 * n, 0.0));
    for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) {
        dist[i][j] = ali->dist(i, j);
        dist[j][i] = dist[i][j];
    }
    vector<vector<double> > q(2 * n, vector<double>(2 * n, 0.0));
    vector<double> r(2 * n, 0.0);

    int m = 0;
    unsigned int i,j,k,f,g;
    while (true) {
        for (i = 0; i < to_merge.size(); i++) {
            r[to_merge[i]] = 0;
            for (j = 0; j < to_merge.size(); j++) r[to_merge[i]] += dist[to_merge[i]][to_merge[j]];
        }
        if (to_merge.size() == 3) break;
        for (i = 0; i < to_merge.size(); i++)
            for (j = 0; j < i; j++)
                q[to_merge[i]][to_merge[j]] = (n - m - 2) * dist[to_merge[i]][to_merge[j]] - r[to_merge[i]] - r[to_merge[j]];
        f = 1;
        g = 0;
        for (i = 0; i < to_merge.size(); i++)
            for (j = 0; j < i; j++)
                if (q[to_merge[i]][to_merge[j]] < q[to_merge[f]][to_merge[g]]) {
                    f = i;
                    g = j;
                }
        Node *node = nodes[n+m];
        node->set_left(nodes[to_merge[f]]);
        node->set_right(nodes[to_merge[g]]);
        nodes[to_merge[f]]->set_father(node);
        nodes[to_merge[g]]->set_father(node);
        nodes[to_merge[f]]->set_dist_father(
            0.5 * dist[to_merge[f]][to_merge[g]] + 0.5 / (n - m - 2) * (r[to_merge[f]] - r[to_merge[g]]));
        nodes[to_merge[g]]->set_dist_father(dist[to_merge[f]][to_merge[g]] - nodes[to_merge[f]]->get_dist_father());
        for (k = 0; k < to_merge.size(); k++) {
            dist[n + m][to_merge[k]] = 0.5 * (dist[to_merge[f]][to_merge[k]] + dist[to_merge[g]][to_merge[k]] - dist[
                                                  to_merge[f]][to_merge[g]]);
            dist[to_merge[k]][n + m] = dist[n + m][to_merge[k]];
        }
        to_merge[f] = n + m;
        to_merge[g] = to_merge.back();
        to_merge.pop_back();
        m++;
    }
    from_root = vector<Node *>();
    nodes[to_merge[0]]->set_dist_father(
        0.5 * dist[to_merge[0]][to_merge[1]] + 0.5 / (n - m - 2) * (r[to_merge[0]] - r[to_merge[1]]));
    from_root.push_back(nodes[to_merge[0]]);
    nodes[to_merge[1]]->set_dist_father(dist[to_merge[0]][to_merge[1]] - nodes[to_merge[0]]->get_dist_father());
    from_root.push_back(nodes[to_merge[1]]);
    nodes[to_merge[2]]->set_dist_father(
        0.5 * (dist[to_merge[0]][to_merge[2]] + dist[to_merge[1]][to_merge[2]] - dist[to_merge[0]][to_merge[1]]));
    from_root.push_back(nodes[to_merge[2]]);
    for (i=0;i<3;i++) from_root[i]->set_father(nullptr);
}


