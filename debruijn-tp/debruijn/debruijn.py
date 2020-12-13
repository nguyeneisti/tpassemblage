#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
#
# import argparse
# import os
# import sys
# import networkx as nx
# import matplotlib
# from operator import itemgetter
# import random
# random.seed(9001)
# from random import randint
# import statistics

import argparse
import os
import sys
import statistics
import random
import operator
import networkx as nx
import matplotlib

__author__ = "Delphine NGUYEN"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Delphine NGUYEN"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Delphine NGUYEN"
__email__ = "nguyendelp@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """ Retourne un générateur de séquences"""
    with open(fastq_file,"rt") as monfic:
        for line in monfic:
            print(line)
            yield next(monfic).strip("\n")
            next(monfic)
            next(monfic)


def cut_kmer(read, kmer_size):
    """Retourne un générateur de k-mer"""
    for i in range(0, len(read)+1-kmer_size):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """ Retourne un dictionnaire : clé = k-mer, valeur = nb occurrences de ce k-mer"""
    dict_k = dict()
    gen_sequence = read_fastq(fastq_file)
    kmer_list = []
    for seq in gen_sequence:
        gen_kmer = cut_kmer(seq,kmer_size)
        kmer_list.append(gen_kmer)
    for kmer in kmer_list:
        for k in kmer:
            if k in dict_k:
                dict_k[k] +=1
            else:
                dict_k[k] = 1
    return dict_k

def build_graph(kmer_dict):
    """ Créer un arbre de k-mers préfixes et suffixes """
    drj_graph = nx.DiGraph() #Directed graph ! Donc on utilise DiGraph, pas juste Graph !
    for elm in kmer_dict.keys():
        drj_graph.add_edge(elm[:-1],elm[1:], weight = kmer_dict[elm])
    return drj_graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Retourne un graphe nettoyé des chemins indésirables"""
    for node in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(node)
        elif delete_entry_node:
            graph.remove_nodes_from(node[:-1])
        elif delete_sink_node:
            graph.remove_nodes_from(node[1:])
        else:
            graph.remove_nodes_from(node[1:-1])
    return graph


def std(data):
    """ Retourne l'écart type de la liste donnée en entrée """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """ Retourne un graphe nettoyé des chemins indésirables """
    best_weights_index = []
    best_lengths = []
    best_paths_index = []
    random.seed(9001)
    for weig in range(len(weight_avg_list)):
        if weight_avg_list[weig] == max(weight_avg_list):
            best_weights_index.append(weig)
    for leng in range(len(path_length)):
        #si la longueur a pour indice le même que celui d'un des meilleurs poids
        if leng in best_weights_index:
            best_lengths.append(path_length[leng])
    for index in best_weights_index:
        if path_length[index]==max(best_lengths):
            best_paths_index.append(index)
    best_path_index = random.choice(best_paths_index)
    path_list.pop(best_path_index)
    path_list2 = path_list
    graph = remove_paths(graph,path_list2,delete_entry_node,delete_sink_node)
    return graph

def path_average_weight(graph, path):
    """ Retourne une moyenne des poids d'un chemin d'un graphe"""
    avg = 0
    for i in range(len(path)-1):
        avg += graph[path[i]][path[i+1]]['weight']
    if avg!=0:
        avg = avg/(len(path)-1)
    else:
        return 0
    return avg

def solve_bubble(graph, ancestor_node, descendant_node):
    """ Retourne un graphe nettoyé de la bulle entre le noeud ancetre et le noeud descendant """
    all_paths = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    weight_avg_list = []
    path_length = []
    for p in all_paths:
        path_length.append(len(p))
        weight_avg_list.append(path_average_weight(graph,p))
    return select_best_path(graph, all_paths, path_length, weight_avg_list)


def simplify_bubbles(graph):
    """ Retourne un graphe sans bulle """
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    all_successors = []
    all_predecessors = []
    starting_node = []
    ending_node = []
    for s_node in starting_nodes:
        for e_node in sink_nodes:
            starting_node.append(s_node)
            ending_node.append(e_node)
            all_successors.append(graph.successors(starting_node[0]))
            all_predecessors.append(graph.predecessors(ending_node[0]))
            if len(all_successors)<2:
                for node in all_successors:
                    starting_node.clear()
                    starting_node = list(node)
                    all_successors.clear()
                    all_successors.append(graph.successors(starting_node[0]))
            if len(all_predecessors)<2:
                for node in all_predecessors:
                    ending_node.clear()
                    ending_node = list(node)
                    all_predecessors.clear()
                    all_predecessors.append(graph.predecessors(ending_node[0]))
            if list(nx.all_simple_paths(graph,starting_node[0],ending_node[0])):
                graph = solve_bubble(graph, starting_node[0], ending_node[0])
    return graph


def solve_entry_tips(graph, starting_nodes):
    """Retourne un graphe sans chemin d'entrée indésirable"""
    graph = simplify_bubbles(graph)
    tous_chemins = [] #liste de liste de chemins, les chemins partant des noeuds d'entrée
    tous_chemins2 = []
    path = []
    path_list2 = []
    end_nodes = get_sink_nodes(graph)
    all_successors = []
    for node in starting_nodes:
        path = []
        s_node = node
        while s_node not in end_nodes:
            path.append(s_node)
            all_successors.append(graph.successors(s_node))
            for i in range(len(all_successors)):
                for successor in all_successors[i]:
                    s_node = successor
        path.append(end_nodes[0])
        tous_chemins.append(path)

    length_list = [] #liste des longueurs de chaque chemin
    weight_list = [] #liste des poids moyens pour chaque chemin
    for i in range(len(tous_chemins)):
        weig = path_average_weight(graph,tous_chemins[i])
        weight_list.append(weig)
    for path in tous_chemins:
        length_list.append(len(path))
    #pour chaque chemin de tous_chemins, on a un poids et une longueur
    # On va trouver le meilleur chemin
    sum = []
    maxi = []
    for i in range(len(weight_list)):
        sum.append(weight_list[i]+(length_list[i]*2.5))
    maxi.append(max(sum))
    best_path_index = sum.index(maxi[0]) #indexe du meilleur chemin
    chemin_a_enlever = [] #liste des noeuds à enlever
    meilleur_chemin = tous_chemins[best_path_index]
    tous_chemins.pop(best_path_index)
    similarity_index = 0
    for elm in tous_chemins:
        for i in range(len(elm)):
            if elm[i] not in meilleur_chemin:
                similarity_index+=1
                chemin_a_enlever.append(elm[i])
        chemin_a_enlever.append(elm[similarity_index])
        break
    path_list2.append(chemin_a_enlever)
    graph2 = remove_paths(graph,path_list2,delete_entry_node=True,delete_sink_node=False)
    return graph2


def solve_out_tips(graph, ending_nodes):
    """Retourne un graphe sans chemin de sortie indésirable"""
    graph = simplify_bubbles(graph)
    tous_chemins = []
    tous_chemins2 = []
    path = []
    path_list2 = []
    start_nodes = get_starting_nodes(graph)
    all_successors = []
    for node in ending_nodes:
        path = []
        s_node = node
        while s_node not in start_nodes:
            path.append(s_node)
            all_successors.append(graph.predecessors(s_node))
            for i in range(len(all_successors)):
                for successor in all_successors[i]:
                    s_node = successor
        path.append(start_nodes[0])
        tous_chemins.append(path)
    tous_chemins3 = []
    for e in tous_chemins:
        temp = e[::-1]
        tous_chemins3.append(temp)
    tous_chemins = tous_chemins3
    length_list = []
    weight_list = []
    for path in tous_chemins:
        weight_list.append(path_average_weight(graph,path))
        length_list.append(len(path))

    sum = []
    maxi = []
    for i in range(len(weight_list)):
        sum.append(weight_list[i]+(length_list[i]*2.5))
    maxi.append(max(sum))
    best_path_index = sum.index(maxi[0])
    chemin_a_enlever = []
    meilleur_chemin = tous_chemins[best_path_index]
    tous_chemins.pop(best_path_index)
    similarity_index = 0
    for elm in tous_chemins:
        for i in range(len(elm)):
            if elm[i] not in meilleur_chemin:
                similarity_index+=1
                chemin_a_enlever.append(elm[i])
        chemin_a_enlever.append(elm[similarity_index])
        break
    path_list2.append(chemin_a_enlever)
    graph2 = remove_paths(graph,path_list2,delete_entry_node=True,delete_sink_node=False)
    return graph2

def get_starting_nodes(graph):
    """ Retourne une liste de noeuds d'entrée """
    nodes_list = []
    for node in graph.nodes():
        if not graph.pred[node]:
            nodes_list.append(node)
    return nodes_list

def get_sink_nodes(graph):
    """ Retourne une liste de noeuds de sortie """
    nodes_list = []
    for node in graph.nodes():
        if not graph.succ[node]:
            nodes_list.append(node)
    return nodes_list

def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Retourne une liste de tuples (contig, taille du contig).
    Un contig est un chemin entre des noeuds d'entrée et des noeuds de sortie.
    """
    contigs_list = [] # liste de tuples : [(c1, len(c1)), (c2, len(c2)), ...]
    for s_node in starting_nodes:
        result = []
        for e_node in ending_nodes:
            paths_list = () #tuple (contig, taille du contig)
            result = []
            paths = nx.all_simple_paths(graph, source=s_node, target=e_node)
            for p in paths:
                for elm in p:
                    result.append(elm[0])
                result.append(elm[1])
            paths_list = (''.join(result), len(result))
            contigs_list.append(paths_list)
    return contigs_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    """ Sauvegarder les contigs dans un fichier texte"""
    file = open(output_file,"w+")
    for i in range(len(contigs_list)):
        fasta = fill(contigs_list[i][0],80)
        file.write(">contig_%d len=%d\n" %(i,contigs_list[i][1]))
        file.write(fasta)
        file.write("\n")
    file.close()

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    fastq_file = args.fastq_file
    k_mer_size = args.kmer_size
    output_file = args.output_file
    dictio = build_kmer_dict(fastq_file,k_mer_size)
    graph = build_graph(dictio)
    graph2 = simplify_bubbles(graph)
    starting_nodes = get_starting_nodes(graph2)
    sink_nodes = get_sink_nodes(graph2)
    graph3 = solve_entry_tips(graph2,starting_nodes)
    graph3 = solve_out_tips(graph3,sink_nodes)
    contigs_list = get_contigs(graph3, starting_nodes, sink_nodes)
    save_contigs(contigs_list, output_file)
    draw_graph(graph3,"graph.png")

if __name__ == '__main__':
    main()
