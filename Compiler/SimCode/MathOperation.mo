/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2016, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

encapsulated package MathOperation
" file:        MathOperation.mo
  package:     MathOperation
  description: create math operation list to use for a adolc trace

  The entry points to this module are the function createOperationData
"


// public imports
public
import Absyn;
import BackendDAE;
import DAE;
import Values;
import SimCode;
import SimCodeVar;

// protected imports
protected
import BackendDump;
import Config;
import ComponentReference;
import DAEUtil;
import Debug;
import Error;
import ErrorExt;
import ExecStat;
import Expression;
import ExpressionDump;
import Flags;
import TaskSystemDump;
import SimCodeUtil;
import System;
import Util;




// data types
public uniontype MathOperator
  record ASSIGN_ACTIVE
  end ASSIGN_ACTIVE;
  record ASSIGN_PARAM
  end ASSIGN_PARAM;
  record ASSIGN_PASSIVE
  end ASSIGN_PASSIVE;
  record PLUS
    Boolean isActive;
  end PLUS;
  record MINUS
    Boolean isActive;
  end MINUS;
  record MUL
    Boolean isActive;
  end MUL;
  record DIV
    Boolean isActive;
  end DIV;
  record POW
  end POW;
  record UNARY_NEG
  end UNARY_NEG;
  record UNARY_CALL
    Absyn.Ident ident;
  end UNARY_CALL;
end MathOperator;

public uniontype Operand
  record OPERAND_VAR
    SimCodeVar.SimVar variable;
  end OPERAND_VAR;
  record OPERAND_CONST
    DAE.Exp const;
  end OPERAND_CONST;
end Operand;

public uniontype Operation
  record OPERATION
    list<Operand> operands;
    MathOperator operator;
    Operand result;
  end OPERATION;
end Operation;

public uniontype OperationData
  record OPERATIONDATA
    list<Operation> operations;
  end OPERATIONDATA;
end OperationData;


/* main entry point for that module */
public function createOperationData
  input list<SimCode.SimEqSystem> inEquations;
  input SimCode.HashTableCrefToSimVar crefToSimVarHT;
  output Option<OperationData> outOperationData;
protected
  list<Operation> operations;
  constant Boolean debug = true;
algorithm
  try
    operations := {};
    if debug then
      print("createOperationData equations input: \n");
      print(Tpl.tplString3(TaskSystemDump.dumpEqs, inEquations, 0, false));
    end if;
    for eq in inEquations loop
      () := matchcontinue eq
      local
        Integer index;
        DAE.Exp exp;
        DAE.ComponentRef cref;
        list<Operation> rest;
        Operation op;
        list<Operand> operands;
        Operand assignOperand, simVarOperand;
        SimCodeVar.SimVar simVar;
        case SimCode.SES_SIMPLE_ASSIGN(index = index, exp = exp, cref = cref) equation
          //operands = {};
          simVar = BaseHashTable.get(cref, crefToSimVarHT);
          simVarOperand = OPERAND_VAR(simVar);
          (assignOperand, operations) = collectOperationsForExp(exp, crefToSimVarHT, operations, 0);
          //print("Done with collectOperationsForExp\n");
          if isTmpOperand(assignOperand) then
            op::rest = operations;
            op = replaceOperationResult(op, simVarOperand);
            operations = op::rest;
          else
            op = OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
            operations = op::operations;
          end if;


        then ();
        //else ();
      end matchcontinue;
    end for;
    operations := listReverse(operations);
    outOperationData := SOME(OPERATIONDATA(operations));
  else
    outOperationData := NONE();
  end try;
end createOperationData;

protected function collectOperationsForExp
  input DAE.Exp inExp;
  input SimCode.HashTableCrefToSimVar crefToSimVarHT;
  input list<Operation> inOps;
  input Integer tmpOffset;
  output Operand lastResult;
  output list<Operation> outOps;
algorithm
  (_, ({lastResult}, outOps, _, _)) := Expression.traverseExpBottomUp(inExp,collectOperation,({}, inOps, tmpOffset, crefToSimVarHT));
end collectOperationsForExp;

protected function collectOperation
  input DAE.Exp inExp;
  input tuple< list<Operand>, list<Operation>, Integer, SimCode.HashTableCrefToSimVar> inTpl;
  output DAE.Exp outExp;
  output tuple<list<Operand>, list<Operation>, Integer, SimCode.HashTableCrefToSimVar> outTpl;
algorithm
  (outExp, outTpl) := matchcontinue (inExp, inTpl)
    local
      DAE.ComponentRef cref;
      SimCode.HashTableCrefToSimVar crefToSimVarHT;
      SimCodeVar.SimVar resVar;
      list<Operand> opds, rest;
      list<Operation> ops;
      Operation operation;
      Operand result;
      DAE.Exp e1,e2;
      DAE.Type ty;
      DAE.Operator op;
      Integer tmpIndex;
      Operand opd1, opd2;
      Absyn.Ident ident;
      String str;

    case (e1 as DAE.RCONST(), (opds, ops, tmpIndex, crefToSimVarHT)) equation
      opds = OPERAND_CONST(e1)::opds;
    then
      (inExp, (opds, ops, tmpIndex, crefToSimVarHT));

    case (DAE.CREF(componentRef=cref), (opds, ops, tmpIndex, crefToSimVarHT)) equation
      resVar = BaseHashTable.get(cref, crefToSimVarHT);
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, tmpIndex, crefToSimVarHT));

    case (DAE.BINARY(operator = op), (opds, ops, tmpIndex, crefToSimVarHT)) equation
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, tmpIndex);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, tmpIndex, crefToSimVarHT));

    case (DAE.CALL(path=Absyn.IDENT("DIVISION"), attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, tmpIndex, crefToSimVarHT)) equation
      op = DAE.DIV(ty);
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, tmpIndex);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, tmpIndex, crefToSimVarHT));

    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, tmpIndex, crefToSimVarHT))
      guard isMathFunction(ident)
    equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, tmpIndex, crefToSimVarHT));

    else
      (inExp, inTpl);
  end matchcontinue;
end collectOperation;

protected function createBinaryOperation
  input DAE.Operator operator;
  input list<Operand> inOpds;
  input Integer inIndex;
  output Operation op;
  output Operand result;
  output Integer outIndex;
algorithm
  _ := match operator
  local
    SimCodeVar.SimVar resVar;
    DAE.Type ty;
    list<Operation> extraOps;
    Boolean isActive, isCommuted;
    Operand opd1, opd2;
    list<Operand> opds;
    DAE.Exp exp;
    case DAE.ADD(ty) equation
      (resVar, outIndex) = createSimTmpVar(inIndex, ty);
      result = OPERAND_VAR(resVar);
      (opds, extraOps, isActive, isCommuted) = checkOperand(inOpds);
      op = OPERATION(opds, PLUS(isActive), result);
    then ();
    case DAE.SUB(ty) equation
      //print("Start SUB\n");
      (resVar, outIndex) = createSimTmpVar(inIndex, ty);
      result = OPERAND_VAR(resVar);
      (opds, extraOps, isActive, isCommuted) = checkOperand(inOpds);
      //print("done with checkOperands isCommuted: " + boolString(isCommuted) + " isActive: " + boolString(isActive) + "\n");
      if not isCommuted then
        op = OPERATION(opds, MINUS(isActive), result);
      else
        (opd1 as OPERAND_CONST(exp))::opd2::{} = opds;
        op = OPERATION({OPERAND_CONST(Expression.negate(exp)), opd2}, PLUS(isActive), result);
      end if;
    then ();
    case DAE.MUL(ty) equation
      (resVar, outIndex) = createSimTmpVar(inIndex, ty);
      result = OPERAND_VAR(resVar);
      (opds, extraOps, isActive, isCommuted) = checkOperand(inOpds);
      op = OPERATION(opds, MUL(isActive), result);
    then ();
    case DAE.DIV(ty) equation
      //print("Start DIV\n");
      (resVar, outIndex) = createSimTmpVar(inIndex, ty);
      result = OPERAND_VAR(resVar);
      (opds, extraOps, isActive, isCommuted) = checkOperand(inOpds);
      //print("done with checkOperands\n");
      if not isCommuted then
        op = OPERATION(opds, DIV(isActive), result);
      else
        (opd1 as OPERAND_CONST(exp))::opd2::{} = opds;
        op = OPERATION({OPERAND_CONST(Expression.invertReal(exp)), opd2}, MUL(isActive), result);
      end if;
    then ();
    else
      fail();
  end match;
end createBinaryOperation;

protected function createSimTmpVar
 input Integer inIndex;
 input DAE.Type inType;
 output SimCodeVar.SimVar simvar;
 output Integer nextInt;
protected
  DAE.ComponentRef cref;
algorithm
  cref := ComponentReference.makeCrefIdent("$tmpOpVar"+"_"+intString(inIndex), inType, {});
  simvar := SimCodeUtil.makeTmpRealSimCodeVar(cref, BackendDAE.TMP_SIMVAR(), inIndex);
  nextInt := inIndex + 1;
end createSimTmpVar;

protected function checkOperand
  input list<Operand> inOpds;
  output list<Operand> outOpds;
  output list<Operation> extraOps = {};
  output Boolean isActive;
  output Boolean isCommuted;
protected
  Operand opd1, opd2;
algorithm
  opd1::opd2::{} := inOpds;
  _ := match(opd1, opd2)

    case (OPERAND_CONST(_), OPERAND_VAR(_)) equation
      outOpds = inOpds;
      isActive = false;
      isCommuted = false;
    then ();

    case (OPERAND_VAR(_), OPERAND_CONST(_)) equation
      outOpds = {opd2,opd1};
      isActive = false;
      isCommuted = true;
    then ();

    case (OPERAND_VAR(_), OPERAND_VAR(_)) equation
      outOpds = inOpds;
      isActive = true;
      isCommuted = false;
    then ();

  end match;
end checkOperand;

protected function isTmpOperand
  input Operand inOperand;
  output Boolean outBoolean;
algorithm
  outBoolean := matchcontinue inOperand
  local
    SimCodeVar.SimVar simVar;
    case OPERAND_VAR(simVar as SimCodeVar.SIMVAR()) equation
      BackendDAE.TMP_SIMVAR() = simVar.varKind;
    then true;

    else false;
  end matchcontinue;
end isTmpOperand;

protected function replaceOperationResult
  input Operation inOp;
  input Operand inOpd;
  output Operation outOp = inOp;
algorithm
  try
    outOp.result := inOpd;
  else
    print(" function replaceOperationResult failed\n");
  end try;
end replaceOperationResult;


protected function isMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"sin", "cos", "exp", "log", "sqrt", "asin", "acos", "asinh", "acosh", "atanh"};
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inName) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isMathFunction;

///////////////////////
/* Dumping functions */
///////////////////////
public function dumpOperationData
  input Option<OperationData> operationData;
protected
  OperationData opData;
  Integer i = 0;
algorithm
  if isSome(operationData) then
    opData := Util.getOption(operationData);
    print("### OpeartionData ###\n");
    for op in opData.operations loop
      i := i+1;
      print(intString(i) + " " + printOperationStr(op) + "\n");
    end for;

  else
    print("No operationData.\n");
  end if;
end dumpOperationData;

protected function printOperationStr
  input Operation inOp;
  output String outString;
algorithm
  outString := "result: " + printOperandStr(inOp.result) + ": ";
  outString := outString + printOperatorStr(inOp.operator) + " ? ";
  for operand in inOp.operands loop
    outString := outString + printOperandStr(operand) + " ";
  end for;
end printOperationStr;

protected function printOperandStr
  input Operand inOp;
  output String outString;
algorithm
  outString := match inOp
  local
    SimCodeVar.SimVar simVar;
    DAE.Exp exp;
    DAE.ComponentRef name;
    Integer index;

    case OPERAND_VAR(simVar as SimCodeVar.SIMVAR(name = name, index = index))
    then "VAR(" + ComponentReference.printComponentRefStr(name) + " " + intString(index) + ") ";

    case OPERAND_CONST(exp)
    then "CONST(" + ExpressionDump.printExpStr(exp) + ") ";

  end match;
end printOperandStr;

protected function printOperatorStr
  input MathOperator inOp;
  output String outString;
algorithm
  outString := match inOp
  local
    Absyn.Ident ident;
    Boolean b;
    case ASSIGN_ACTIVE()
    then "ASSIGN_ACTIVE";

    case ASSIGN_PARAM()
    then "ASSIGN_PARAM";

    case ASSIGN_PASSIVE()
    then "ASSIGN_PASSIVE";

    case PLUS(b)
    then "PLUS" + (if b then "_ACTIVE" else "_PASSIVE");

    case MINUS(b)
    then "MINUS" + (if b then "_ACTIVE" else "_PASSIVE");

    case MUL(b)
    then "MUL" + (if b then "_ACTIVE" else "_PASSIVE");

    case DIV(b)
    then "DIV" + (if b then "_ACTIVE" else "_PASSIVE");

    case POW()
    then "POW";

    case UNARY_NEG()
    then "UNARY_NEG";

    case UNARY_CALL(ident)
    then "UNARY_CALL_" + ident;

  end match;
end printOperatorStr;


annotation(__OpenModelica_Interface="backend");
end MathOperation;
