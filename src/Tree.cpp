//
// Created by Didelot, Xavier on 29/09/2025.

#include "Tree.h"
#include "Alignment.h"

#include <iostream>

/**
 *
 * @param n
 */
Tree::Tree(int n): n(n), nodes(vector<Node *>(n+n-1)) {
    for (int i = 0; i < n; i++) {
        Node *node = new Node();
        node->set_name("seq" + std::to_string(i + 1));
        node->set_dist_father(1.0);
        nodes[i]=node;
    }

    for (int i = 0; i < n-1; i++) {
        Node *node = new Node();
        nodes[n+i]=node;
        node->set_left(nodes[i]);
        nodes[i]->set_father(node);
        node->set_dist_father(1.0);
        if (i>0) {
            node->set_father(nodes[n+i-1]);
            nodes[n+i-1]->set_right(node);
        }
    }
    nodes[n+n-2]->set_right(nodes[n-1]);
    nodes[n-1]->set_father(nodes[n+n-2]);
    nodes[n-1]->set_dist_father(1.0);
}

string Tree::to_string() {
    return nodes[n]->to_string() + ";";
}

Tree::~Tree() {
    for (unsigned int i = 0; i < nodes.size(); i++) delete nodes[i];
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
    if (n * 2 - 1 - nodes.size() != 0) return false;
    for (unsigned int i = 0; i < nodes.size(); i++) {
        if (nodes[i]->get_left() != nullptr && nodes[i]->get_left()->get_father() != nodes[i]) return false;
        if (nodes[i]->get_right() != nullptr && nodes[i]->get_right()->get_father() != nodes[i]) return false;
    }
    return true;
}

void Tree::add_mutations(int l, double theta) {
    delete nodes[n]->get_seq();
    nodes[n]->set_seq(new vector<uint8_t>(l));
    uniform_int_distribution unif4(0, 3);
    for (int i = 0; i < l; i++) {
        nodes[n]->get_seq()->at(i) = Alignment::nuc_from_int(unif4(rng));
    }
    nodes[n]->get_left ()->add_mutations(nodes[n]->get_seq(), theta);
    nodes[n]->get_right()->add_mutations(nodes[n]->get_seq(), theta);
}

void Tree::randomise(double avg_bra_len) {
    vector<Node *> to_merge = vector<Node *>(n);
    exponential_distribution rexp(avg_bra_len);
    for (int i = 0; i < n; i++) {
        to_merge[i]=nodes[i];
        to_merge[i]->set_dist_father(rexp(rng));
    }
    int m=n+1;
    while (to_merge.size() > 2) {
        uniform_int_distribution unif(0, static_cast<int>(to_merge.size() - 1));
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
    nodes[n]->set_left (to_merge[0]);
    nodes[n]->set_right(to_merge[1]);
    to_merge[0]->set_father(nodes[n]);
    to_merge[1]->set_father(nodes[n]);
}

void Tree::nj(Alignment *ali) {
    vector<int> to_merge = vector<int>(n);
    for (int i = 0; i < n; i++) to_merge[i]=i;
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
    nodes[to_merge[0]]->set_dist_father(
        0.5 * dist[to_merge[0]][to_merge[1]] + 0.5 / (n - m - 2) * (r[to_merge[0]] - r[to_merge[1]]));
    nodes[to_merge[1]]->set_dist_father(dist[to_merge[0]][to_merge[1]] - nodes[to_merge[0]]->get_dist_father());
    nodes[to_merge[2]]->set_dist_father(
        0.5 * (dist[to_merge[0]][to_merge[2]] + dist[to_merge[1]][to_merge[2]] - dist[to_merge[0]][to_merge[1]]));
    nodes[n+m]->set_left (nodes[to_merge[0]]);
    nodes[n+m]->set_right(nodes[to_merge[1]]);
    nodes[to_merge[0]]->set_father(nodes[n+m]);
    nodes[to_merge[1]]->set_father(nodes[n+m]);
    nodes[n+m+1]->set_left (nodes[to_merge[2]]);
    nodes[n+m+1]->set_right(nodes[n+m]);
    nodes[to_merge[2]]->set_father(nodes[n+m+1]);
    nodes[n+m]->set_father(nodes[n+m+1]);
    nodes[n+m]->set_dist_father(0.0);
    swap(nodes[n+m+1],nodes[n]);
}

void Tree::set_seqs(Alignment *ali) {
    for (int i=0;i<ali->get_n();i++) get_node(i)->set_seq(ali->get_seq(i));
}

void Tree::ancestral_reconstruction() {
    unsigned int l=nodes[0]->get_seq()->size();
    for (unsigned int i=n;i<nodes.size();i++) {
        delete nodes[i]->get_seq();
        nodes[i]->set_seq(new vector<uint8_t>(l,0));
    }
    for (unsigned int j=0;j<l;j++) {
        uint8_t combined=0u;
        for (int i=0;i<n;i++) {
            if (Alignment::nuc_is_missing(nodes[i]->get_seq()->at(j))) continue;
            combined=combined|nodes[i]->get_seq()->at(j);
        }
        if (Alignment::nuc_is_single(combined)) {
            for (int i=0;i<n;i++) Alignment::nuc_set_as(&nodes[i]->get_seq()->at(j),combined);
            for (unsigned int i=n;i<nodes.size();i++) nodes[i]->get_seq()->at(j)=combined;
            continue;
        }
        //Fitch algo to go here since site is polymorphic
    }
}


