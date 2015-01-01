#!/usr/bin/env python
# encoding=utf8
#
# Created on:  Mar 13, 2013
#

"""A model for scientific collaboration networks considering tie strength.
For detailed description, please see http://arxiv.org/abs/1401.5027.
"""

__author__ = 'Qing Ke'
__email__ = 'keqing.echolife@gmail.com'


import networkx as nx
import os
import random

from collections import Counter
from itertools import combinations
from weighted_sample import random_weighted_sample_no_replacement, weighted_choice


class Student:
    """
    """
    def __init__(self, author_id, init_exp=1):
        """
        Parameters
        ----------
        author_id : integer

        init_exp : integer
            Initial expertise for the student
        """
        self.author_id = author_id
        self.exp = init_exp

    def inc_experience(self):
        """Increase the expertise by 1.
        """
        self.exp += 1


class Advisor:
    """
    """
    def __init__(self, author_id):
        """
        Parameters
        ----------
        author_id : integer
        """
        self.author_id = author_id


class Group:
    """
    """
    def __init__(self, advisor):
        """
        Parameters
        ----------
        advisor : Advisor
            The advisor for the group
        """
        self.advisor = advisor
        self.students = [] 

    def add_student(self, student):
        """Add a student to the group.

        Parameters
        ----------
        student : Student
        """
        self.students.append(student)


class CA:
    """
    """
    def __init__(self, new_stu_num, grp_ppr_num, thr, inactive_p, ratio):
        """Initiate the scientific collaboration system.

        Parameters
        ----------
        new_stu_num : integer
            Number of new students joining each group at each time step

        grp_ppr_num : integer
            Number of intra-group papers at each time step

        thr : integer
            Expertise threshold for students to graduate

        inactive_p : float
            Inactive probability for students after graduation

        ratio : integer
            Ratio of inter- and intra-group papers
        """
        self.groups = []  # all the groups
        self.b_out = open('model_bipartite.txt', 'w') # author-paper bipartite network
        self.network = nx.Graph()  # collaboration network
        self.total_id = 0  # author (advisor and students) id
        self.new_stu_num = new_stu_num
        self.grp_ppr_num = grp_ppr_num
        self.thr = thr
        self.inactive_p = inactive_p
        self.ratio = ratio
        self.colla_dict = {}
        self.colla_total = set()
        self.start(1, 2)

    def get_next_id(self):
        """Return the next available author id.
        """
        self.total_id += 1
        self.network.add_node(self.total_id)
        return self.total_id

    def new_group(self, advisor_id):
        """Create a new group.

        Parameters
        ----------
        advisor_id : integer
            Id for the group advisor
        """
        advisor = Advisor(advisor_id)
        group = Group(advisor)
        self.groups.append(group)
        return len(self.groups) - 1

    def start(self, group_num, group_size):
        """Initial configuration.

        Parameters
        ----------
        group_num : integer
            Number of initial groups

        group_size : integer
            Number of initial members in each group including 1 advisor
        """
        for i in range(group_num):  # each group
            grp_index = self.new_group( self.get_next_id() )
            assert grp_index == i
            for j in range(group_size - 1):  # with group_size-1 students
                student = Student( self.get_next_id() )
                self.groups[i].add_student(student)

    def evolve(self, time_step):
        """Evolve the system for a certain time steps.

        Parameters
        ----------
        time_step : integer
            Number of time steps
        """
        num_group, num_stud = [], []
        for i in range(time_step):
            print 'Time step', i
            self.each_step()
            num_group.append( len(self.groups) )
            total = sum( len(grp.students) for grp in self.groups )
            num_stud.append(total)
        return num_group, num_stud

    def each_step(self):
        """At each time step.
        """
        # intra-group
        self.publish_grp_paper()

        # inter-group
        if len(self.groups) > 1:
            for j in range(self.grp_ppr_num * self.ratio):
                for i in range( len(self.groups) ):
                    if random.random() < 0.4:
                        self.collaborate_fixed(i)

        # new group
        new_groups_advisor_id = self.check_exp()
        for advisor_id in new_groups_advisor_id:
            self.new_group(advisor_id)

        # new students
        self.join_new_stu()

    def publish_grp_paper(self):
        """Intra-group paper. One author is the group advisor,
        the rest are chosen from students based on their expertise.
        """
        for group in self.groups: # each group 
            for i in range(self.grp_ppr_num):  # publishes grp_ppr_num paper
                # publishing probability controls the ratio of total number of authors to that of papers
                if random.random() < 0.4:
                    authors = [group.advisor.author_id]
                    self.choose_students(group.students, authors) # choose co-authors, in place
                    if len(authors) == 1: # at least 2 authors, just in case
                        continue
                    self.b_out.write( '\t'.join([str(a) for a in authors]) + '\n' )
                    for pair in combinations(authors, 2):
                        self.cumu_weight( pair, len(authors) )

    def choose_students(self, students, authors):
        """Choose co-authors randomly from students list based on their expertise, 
        and append to authors list. If the number of students are not greater than
        l-1, all the students are chosen.

        Note that the elements in authors are author_id.

        Parameters
        ----------
        students : list
            List of students

        authors : list
            List of authors
        """
        # The number of authors in each paper is obtained from the APS data
        choices = [32,28,16,9,7,4,3,2,1]
        idx_to_thr = {0:2, 1:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10}
        idx = weighted_choice(choices)
        ath_num = idx_to_thr[idx]

        if len(students) < ath_num - 1: # not enough co-authors
            for stu in students:
                authors.append(stu.author_id)
        else:
            items = []
            for stu in students:
                items.append( (stu.exp, stu.author_id) )
            co_aut = list( random_weighted_sample_no_replacement(items, ath_num-1) )
            authors.extend(co_aut)

    def collaborate_fixed(self, grp_one_index):
        """Inter-group papers. Two authors are the two groups advisors and others randomly
        select from the pool of students of the two groups based on expertise.

        Parameters
        ----------
        grp_one_index : integer
            Index for the group
        """
        grp_one = self.groups[grp_one_index]
        grp_one_leader = grp_one.advisor.author_id

        # choose group to collaborate
        if grp_one_index in self.colla_dict: # established collaboration
            grp_two_index = self.colla_dict[grp_one_index]
        else: # no collaborative group
            total = set( range(len(self.groups)) )
            choices = list(total - self.colla_total)
            choices.remove(grp_one_index)
            if len(choices) == 0:
                return
            grp_two_index = random.choice(choices)
            self.colla_dict[grp_one_index] = grp_two_index
            self.colla_dict[grp_two_index] = grp_one_index
            self.colla_total.add(grp_one_index)
            self.colla_total.add(grp_two_index)

        # collaborate
        grp_two = self.groups[grp_two_index]
        grp_two_leader = grp_two.advisor.author_id
        authors = [grp_one_leader, grp_two_leader]

        items = []
        for stu in grp_one.students:
            items.append( (stu.exp, stu.author_id) )
        for stu in grp_two.students:
             items.append( (stu.exp, stu.author_id) )

        # based on APS data
        choices = [32,28,16,9,7,4,3,2,1]
        idx_to_thr = {0:2, 1:3, 2:4, 3:5, 4:6, 5:7, 6:8, 7:9, 8:10}
        idx = weighted_choice(choices)
        ath_num = idx_to_thr[idx]
        co_aut = list( random_weighted_sample_no_replacement(items, ath_num-2) )
        authors.extend(co_aut)

        self.b_out.write( '\t'.join([str(a) for a in authors]) + '\n' )
        for pair in combinations(authors, 2):
            self.cumu_weight( pair, len(authors) )

    def check_exp(self):
        """Students whose expertise reach the threshold can be inactive with probability
        inactive_p or lead a new group with probability 1-inactive_p.

        Returns
        -------
        The list of advisors' id of new groups
        """
        created_groups_advisor_id = []
        for group in self.groups:
            removed_stu = []
            for stu in group.students:
                stu.inc_experience()
                if stu.exp >= self.thr:
                    if random.random() > self.inactive_p:
                        created_groups_advisor_id.append(stu.author_id)
                    removed_stu.append(stu)
            for stu in removed_stu:
                group.students.remove(stu)
        return created_groups_advisor_id

    def join_new_stu(self):
        """Each group has new students.
        """
        for group in self.groups:
            for i in range(self.new_stu_num):
                new_stu = Student( self.get_next_id() )
                group.add_student(new_stu)

    def cumu_weight(self, pair, author_num):
        """Cumulate the weight of given edge.

        Parameters
        ----------
        pair : tuple
            An edge

        author_num : integer
            Number of authors
        """
        src, dst = pair[0], pair[1]
        if dst in self.network[src]: # neighbors
            self.network[src][dst]['weight'] += 1.0/(author_num-1)
        else:
            self.network.add_edge( src, dst, weight=1.0/(author_num-1) )

    def analyze_network(self, out_path):
        """Analyze the collaboration network after the completion of evolving.

        Parameters
        ----------
        out_path : string
            Output file path
        """
        lcc = nx.connected_component_subgraphs(self.network)[0]
        nx.write_weighted_edgelist(lcc, out_path)
        total_stud = sum( len(grp.students) for grp in self.groups )
        return len(self.groups), total_stud, lcc.number_of_nodes(), lcc.number_of_edges(), nx.average_clustering(lcc), nx.degree_assortativity_coefficient(lcc)


def main(B):
    """Repeat simulation B times.
    """
    new_stu_num, grp_ppr_num, thr, inactive_p, ratio = 1, 1, 7, 0.8, 1
    time_step = 65

    f_str = open('structural_results.txt', 'w')
    f_time = open('grp_stud_time_results.txt', 'w')
    for i in range(B):
        print 'Simulation', i
        out_path = 'model_colla_network_ratio_' + str(ratio) + '_' + str(i) + '.wedgelist'
        colla = CA(new_stu_num, grp_ppr_num, thr, inactive_p, ratio)
        num_group, num_stud = colla.evolve(time_step)
        f_time.write( str(num_group) + '\t' + str(num_stud) + '\n' )

        G, STUD, N, M, CC, A = colla.analyze_network(out_path)
        out = [ str(G), str(STUD), str(N), str(M), str(2.0*M/N), str(CC), str(A) ]
        f_str.write( '\t'.join(out) + '\n')


if __name__ == '__main__':
    main(5)
