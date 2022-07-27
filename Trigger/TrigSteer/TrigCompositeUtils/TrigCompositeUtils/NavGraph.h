/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TrigCompositeUtils_NavGraph_h
#define TrigCompositeUtils_NavGraph_h

#include "AthContainers/AuxElement.h"
#include "xAODTrigger/TrigCompositeContainer.h"
#include "xAODTrigger/TrigCompositeAuxContainer.h"


namespace TrigCompositeUtils {

  /**
   * @class NavGraphNode
   * @brief Transient utility class to represent a node in a graph (m_decisionObject), and a vector of edges (m_filteredSeeds)
   * to other nodes which are parents of this node.
   **/
  class NavGraphNode {
    public:

      /**
       * @brief Construct a NavGraphNode shadowing a node in the full xAOD navigation graph
       * @param[in] me The Decision object node from the full xAOD navigation graph which this object is representing.
       **/
      NavGraphNode(const Decision* me);

      /**
       * @brief Form an edge in the graph from this node to another one.
       * @param[in] to The "parent" or "seed" Decision object from the perspective of this Node's shadowed Decision object. Mutable to allow two-way linking.
       * @return True if a new edge was added. False if this was a duplicated call to add this edge.
       **/
      bool linksTo(NavGraphNode* to);

      /**
       * @brief Forget about any graph edges to the supplied node. Forgets both child and seed (a.k.a. parent) linking
       * @param[in] Node to un-link
       **/
      void dropLinks(NavGraphNode* node);

      /**
       * @brief Flag this node as one to keep when the thin() operation is performed
       **/
      void keep();

      /**
       * @brief Reset the keep flag to false upon finishing thinning
       **/
      void resetKeep() ;

      /**
       * @return If the keep flag was set
       **/
      bool getKeep() const;

      /**
       * @brief Return a const pointer to the Decision object node which this NavGraphNode is shadowing.
       **/
      const Decision* node() const;

      /**
       * @brief Return a vector of const pointers to the Decision object nodes which this NavGraphNode seeds from. A.k.a its parents.
       * Note: NavGraph is used to represent a sub-graph of the full navigation graph, hence it is expected that
       * the vector of seeds returned from this function may be smaller than the vector of seeds returned from the
       * shadowed xAOD Decision Object.
       **/
      const std::vector<NavGraphNode*>& seeds() const;

      /**
       * @brief Return a vector of const pointers to the Decision object nodes which are the children of this NavGraphNode.
       * Note: The m_decisionObject does not provide such forward-exploring capability.
       **/
      const std::vector<NavGraphNode*>& children() const;

    private:

      /**
       * @brief Internal helper function. Using a vector to preserve pointer ordering, but want the de-duplication behavior of a vector.
       * @return true if toAdd was added, false if it was already contained in the vector.
       **/
      static bool addIfNotDuplicate(std::vector<NavGraphNode*>& container, NavGraphNode* toAdd);

      const Decision* m_decisionObject; //!< The Decision object node which I shadow
      std::vector<NavGraphNode*> m_filteredSeeds; //!< My seeds (edges in the graph), filtered on per-chain requirements.
      std::vector<NavGraphNode*> m_filteredChildren; //!< Two-way linking information, used when thinning the graph.
      bool m_keepFlag; //!< Keep this node when slimming the NavGraph. Needs to be set explicitly

  };

  /**
   * @class NavGraph
   * @brief Structure to hold a transient Directed Acyclic Graph (DAG) structure.
   * NavGraph is populated from, and forms a sub-graph over the full Run 3 trigger navigation graph in a single event.
   * Requirements on specific chains, and the specification of allowable graph entry-points are considered in the construction of the NavGraph sub-graph.
   * Once one of these sub-graphs is fully populated to a given specification, it is searched by the feature-retrieval code to find features.
   **/
  class NavGraph {

    public:

      /**
       * @brief Construct an empty NavGraph
       **/
      NavGraph();

      /**
       * @brief Add a new NavGraphNode which shadows the xAOD Decision object "node" from the full navigation graph
       * @param[in] node The xAOD Decision object which the new node will shadow. Will not cause duplication if node has already been added.
       * @param[in] ctx The event context.
       * @param[in] comingFrom If not null, used to indicate which xAOD Decision object was the seed of "node". This is used to form an edge in the graph.
       * Alternately, if comingFrom is null then "node" is taken as a final node (one of the locations from which the graph should be explored) and hence is added
       * to the finalNodes vector. 
       **/
      void addNode(const Decision* node, const EventContext& ctx, const Decision* comingFrom = nullptr);

      /**
       * @brief Get all final nodes.
       * @return Vector of final nodes. These are the nodes which were added without any "comingFrom". 
       * To explore the NavGraph fully, one should explore recursively all paths originating from each of the final nodes.
       **/
      std::vector<NavGraphNode*> finalNodes() const;

      /**
       * @brief Get all nodes.
       * @return Vector of all nodes. Including all final, intermediate, and initial nodes.
       **/
      std::vector<NavGraphNode*> allNodes();

      /**
       * @return Total number of nodes in the NavGraph. 
       **/
      size_t nodes() const;

      /**
       * @return Total number of edges in the NavGraph. 
       **/
      size_t edges() const;

      /**
       * @brief Perform thinning. Removing all nodes which are not explicitly flagged as keep(), after having re-wired them out of the graph.
       * @return A vector of the Decision* behind the NavGraphNodes which were thinned from the NavGraph.
       **/
      std::vector<const Decision*> thin();

      /**
       * @brief Helper function. Print the internal graph structure to the terminal.
       * @param[in] log Athena messaging service reference.
       * @param[in] msgLevel Athena messaging service verbosity level.
       **/
      void printAllPaths(MsgStream& log, MSG::Level msgLevel = MSG::VERBOSE) const;

    private:

      /**
       * @brief Take all seeds (parents) of the supplied node and connect them to all the node's children. 
       * Unlink the parents and children from the node.
       * For the case of Parent nodes, P, node to be deleted, N, and Child nodes, C, this function converts from
       *   P      P  P     P      P  P
       *   |      \ /      |      \ /
       *   N   ,   N   ,   N   ,   N
       *   |       |      / \     / \
       *   C       C      C  C    C  C
       * to
       *   P      P  P     P      P   P
       *   |      \ /      |      |\ /|
       *   |   ,   |   ,   |   ,  | | | 
       *   |       |      / \     |/ \|
       *   C       C      C  C    C   C
       * where N is orphaned from the graph, with no parents or children.
       * @param[in] toBeDeleted Node to rewire out of the navigation graph prior to its removal.
       **/
      void rewireNodeForRemoval(NavGraphNode& toBeDeleted);

      /**
       * @bried Internal helper function. Recursively print the graph structure from a single given starting node.
       * @param[in] nav The node to recursively explore.
       * @param[in] level The current depth of recursion. Used to pad output format.
       * @param[in] log Athena messaging service reference.
       * @param[in] msgLevel Athena messaging service verbosity level.
       **/
      void recursivePrintNavPath(const NavGraphNode& nav, size_t level, MsgStream& log, MSG::Level msgLevel) const;

      std::map<const ElementLink<TrigCompositeUtils::DecisionContainer>, NavGraphNode> m_nodes; //!< Map of nodes in the graph. Indexed on the underlying Decision object's ElementLink.
      std::vector<NavGraphNode*> m_finalNodes; //!< Entry points into the navigation graph. When iterating over the graph, start from all of these places.
      size_t m_edges; //!< Statistics on the number of edges, connecting the nodes in the graph.
  };

}

#endif // TrigCompositeUtils_NavGraph_h
