#ifndef __NODEVISITOR_H_
#define __NODEVISITOR_H_

#include <vector>

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/subs.h>
#include <symengine/logic.h>
#include <symengine/matrix.h>

namespace NodeVisitor
{
  using namespace SymEngine;
  using namespace std;

  // This is a visitor class that transverse the expression tree of a SymEngine object, creating a representation
  // of the expression that can be saved and used later to reconstruct the original expression.
  struct Node: public BaseVisitor<Node>
  {
    enum class NodeType
      {
       ADD,
       MUL,
       POW,
       SYM,
       INT,
       FUNC,
      };
    // The type of the SymEngine object represented by this node.
    NodeType type;
    // The name of the SymEngine object represented by this node.
    // Only appliable if its is a SYM, INT or FUNC.
    string name{};
    // The vector of childs
    vector<Node> childs{};

    Node() = default;
    Node(const RCP<const Basic> &b) {b->accept(*this);}
    bool operator==(const Node &n) const
    {
      return type == n.type and name == n.name
        and childs == n.childs;
    }
    bool operator!=(const Node &n) const
    {
      return not (*this == n);
    }

    // Return the expression being represented by the tree whose root is the current Node.
    RCP<const Basic> to_basic() const;


    // This helper method populates the vector 'childs' of the current Node.
    void parse(const Basic & expr)
    {
      for (const auto & e : expr.get_args())
        {
          Node n;
          e->accept(n);
          childs.emplace_back(move(n));
        }
    }

    void bvisit(const Add & expr)
    {
      type = NodeType::ADD;
      parse(expr);
    }

    void bvisit(const Mul & expr)
    {
      type = NodeType::MUL;
      parse(expr);
    }

    void bvisit(const Pow & expr)
    {
      type = NodeType::POW;
      parse(expr);
    }

    void bvisit(const Symbol & expr)
    {
      type = NodeType::SYM;
      name = expr.get_name();
    }

    void bvisit(const Integer & expr)
    {
      type = NodeType::INT;
      name = expr.__str__();
    }

    void bvisit(const FunctionSymbol & expr)
    {
      type = NodeType::FUNC;
      name = expr.get_name();
      parse(expr);
    }

    template<typename T>
    void bvisit(const T &e)
    {
      throw runtime_error{"Can not make a tree expression from " + e.__str__()};
    }

    template<typename S>
    void serialize(S& s){
      s(type, name, childs);
    }
  };

  RCP<const Basic> Node::to_basic() const
  {
    RCP<const Basic> res;
    switch (type)
      {
      case NodeType::ADD:
        res = zero;
        for (const auto & n : childs)
          res = add(res, n.to_basic());
        break;
      case NodeType::MUL:
        res = one;
        for (const auto & n : childs)
          res = mul(res, n.to_basic());
        break;
      case NodeType::POW:
        res = pow(childs[0].to_basic(), childs[1].to_basic());
        break;
      case NodeType::SYM:
        res = symbol(name);
        break;
      case NodeType::INT:
        res = integer(name);
        break;
      case NodeType::FUNC:
        auto fx{function_symbol(name,childs[0].to_basic())};
        res = fx;
        break;

      }
    return res;
  }

  // Given a DenseMatrix, returns "matrix" of Nodes representing each element of the original matrix.
  vector<vector<Node>> matrix2vecnode(const DenseMatrix & dm)
  {
    vector<vector<Node>> vec(dm.nrows());
    for (auto i = 0; i < dm.nrows(); i++)
      for (auto j = 0; j < dm.ncols(); j++)
        vec[i].emplace_back(dm.get(i,j));

    return vec;
  }

  // Given a "matrix" of Nodes, reconstructs the DenseMatrix being represented by these Nodes.
  DenseMatrix vecnode2matrix(const vector<vector<Node>> & vec)
  {
    auto n{vec.size()}, m{vec.at(0).size()};
    DenseMatrix dm{n,m};
    for (auto i = 0; i < n; i++)
      for (auto j = 0; j < m; j++)
        dm.set(i,j, vec.at(i).at(j).to_basic());

    return dm;
  }


}



#endif // __NODEVISITOR_H_
