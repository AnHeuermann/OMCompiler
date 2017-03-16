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

/* TODO:
    - write new entry point for a model
      - define independets and dependents
      - care for functions
      -
    - create operationsDataStmts
    - use workingStateArgs in alls functions
*/


// data types
public uniontype MathOperator
  record ASSIGN_ACTIVE
  end ASSIGN_ACTIVE;
  record ASSIGN_PARAM
  end ASSIGN_PARAM;
  record ASSIGN_PASSIVE
  end ASSIGN_PASSIVE;
  record ASSIGN_IND
  end ASSIGN_IND;
  record ASSIGN_DEP
  end ASSIGN_DEP;
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
  record MODELICA_CALL
    Absyn.Ident ident;
  end MODELICA_CALL;
end MathOperator;

public uniontype Operand
  record OPERAND_VAR
    SimCodeVar.SimVar variable;
  end OPERAND_VAR;
  record OPERAND_CONST
    DAE.Exp const;
  end OPERAND_CONST;
  record OPERAND_TIME
  end OPERAND_TIME;
  record OPERAND_INDEX
    Integer i;
  end OPERAND_INDEX;
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
    Integer maxTmpIndex;
    list<Integer> independents;
    list<Integer> dependents;
  end OPERATIONDATA;
end OperationData;

protected uniontype WorkingStateArgs
  record WORKINGSTATEARGS
    SimCode.HashTableCrefToSimVar crefToSimVarHT;
    list<Absyn.Path> funcNames;
    Integer maxTmpIndex;
    Integer numParameters;
  end WORKINGSTATEARGS;
end WorkingStateArgs;


public function createOperationData
  input list<SimCode.SimEqSystem> inEquations;
  input SimCode.HashTableCrefToSimVar crefToSimVarHT;
  input SimCode.VarInfo varInfo;
  input DAE.FunctionTree functions;
  output Option<OperationData> outOperationData;
protected
  WorkingStateArgs workingArgs;
  Integer numVariables; 
algorithm
  try
    numVariables := 2*varInfo.numStateVars + varInfo.numAlgVars;
    workingArgs := WORKINGSTATEARGS(crefToSimVarHT, {}, numVariables, varInfo.numParams);
    
    // create operation of the equations
    (outOperationData,workingArgs) := createOperationEqns(inEquations, workingArgs);
    // create needed functions
    
  else
    outOperationData := NONE();  
  end try;
end createOperationData;



/* main entry point for that module */
protected function createOperationEqns
  input list<SimCode.SimEqSystem> inEquations;
  input WorkingStateArgs iworkingArgs;
  output Option<OperationData> outOperationData;
  output WorkingStateArgs workingArgs = iworkingArgs;
protected
  list<Operation> operations, tmpOps;
  constant Boolean debug = false;
  Integer tmpIndex;
  list<DAE.Statement> statements;
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
          simVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
          simVarOperand = OPERAND_VAR(simVar);
          (assignOperand, operations, workingArgs) = collectOperationsForExp(exp, operations, workingArgs);
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

        case SimCode.SES_ALGORITHM(index = index, statements=statements) equation
          (tmpOps, workingArgs) = createOperationDataStmts(statements, workingArgs);
          operations = listAppend(tmpOps, operations); 
        then ();
        //else ();
      end matchcontinue;
    end for;
    operations := listReverse(operations);
    outOperationData := SOME(OPERATIONDATA(operations, workingArgs.maxTmpIndex, {}, {}));
  else
    outOperationData := NONE();
  end try;
end createOperationEqns;


public function createOperationDataStmts
  input list<DAE.Statement> inStmts;
  input WorkingStateArgs iworkingArgs;
  output list<Operation> outOperationData;
  output WorkingStateArgs workingArgs = iworkingArgs;
protected
  list<Operation> operations;
  constant Boolean debug = false;
  Integer tmpIndex;
algorithm
  try
    operations := {};
    if debug then
      print("createOperation for statements input: \n");
    end if;
    for smts in inStmts loop
      () := matchcontinue smts
      local
        Integer index;
        DAE.Exp rhs, lhs;
        DAE.ComponentRef cref;
        list<Operation> rest;
        Operation op;
        list<Operand> operands;
        Operand assignOperand, simVarOperand;
        SimCodeVar.SimVar simVar;
        case DAE.STMT_ASSIGN(exp = lhs, exp1 = rhs) equation
          //operands = {};
          cref = Expression.expCref(lhs);
          simVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
          simVarOperand = OPERAND_VAR(simVar);
          (assignOperand, operations, workingArgs) = collectOperationsForExp(rhs, operations, workingArgs);
          //print("Done with collectOperationsForExp\n");
          if isTmpOperand(assignOperand) then
            op::rest = operations;
            op = replaceOperationResult(op, simVarOperand);
            operations = op::rest;
          else
            op = OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
            outOperationData = op::operations;
          end if;
        then ();
        //else ();
      end matchcontinue;
    end for;
  else
    outOperationData := {};
  end try;
end createOperationDataStmts;

protected function collectOperationsForExp
  input DAE.Exp inExp;
  input list<Operation> inOps;
  input WorkingStateArgs iworkingArgs;
  output Operand lastResult;
  output list<Operation> outOps;
  output WorkingStateArgs oworkingArgs;
algorithm
  (_, ({lastResult}, outOps, oworkingArgs)) := Expression.traverseExpBottomUp(inExp,collectOperation,({}, inOps, iworkingArgs));
end collectOperationsForExp;

protected function collectOperation
  input DAE.Exp inExp;
  input tuple< list<Operand>, list<Operation>, WorkingStateArgs> inTpl;
  output DAE.Exp outExp;
  output tuple<list<Operand>, list<Operation>, WorkingStateArgs> outTpl;
algorithm
  (outExp, outTpl) := matchcontinue (inExp, inTpl)
    local
      DAE.ComponentRef cref;
      SimCode.HashTableCrefToSimVar crefToSimVarHT;
      SimCodeVar.SimVar resVar, timeVar, paramVar, extraVar, startVar;
      list<Operand> opds, rest;
      list<Operation> ops;
      Operation operation;
      Operand result;
      DAE.Exp e1,e2;
      list<DAE.Exp> expList;
      DAE.Type ty;
      DAE.Operator op;
      Integer tmpIndex;
      Operand opd1, opd2;
      Absyn.Ident ident;
      String str;
      WorkingStateArgs workingArgs;

    case (e1 as DAE.RCONST(), (opds, ops, workingArgs)) equation
      opds = OPERAND_CONST(e1)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CREF(componentRef=DAE.CREF_IDENT(ident="time"), ty=ty), (opds, ops, workingArgs)) equation
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, ty);
      operation = OPERATION({OPERAND_TIME()}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CREF(componentRef=cref, ty=ty), (opds, ops, workingArgs)) equation
      paramVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      //guard
      BackendDAE.PARAM() = paramVar.varKind;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, ty);
      workingArgs.maxTmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_VAR(paramVar)}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CREF(componentRef=cref), (opds, ops, workingArgs)) equation
      resVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT("der")), (opds, ops, workingArgs)) equation
      OPERAND_VAR(resVar)::rest = opds;
      cref = ComponentReference.crefPrefixDer(resVar.name);
      resVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      opds = OPERAND_VAR(resVar)::rest;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.BINARY(operator = op), (opds, ops, workingArgs)) equation
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, workingArgs.maxTmpIndex);
      workingArgs.maxTmpIndex = tmpIndex;
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT("DIVISION"), attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, workingArgs)) equation
      op = DAE.DIV(ty);
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, workingArgs.maxTmpIndex);
      workingArgs.maxTmpIndex = tmpIndex;
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
      guard isMathFunction(ident)
    equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, ty);
      workingArgs.maxTmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
      guard isTrigMathFunction(ident)
    equation
      opd1::rest = opds;
      (extraVar, tmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, ty);
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      workingArgs.maxTmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1,OPERAND_VAR(extraVar)}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT("start"),  expLst={DAE.CREF(componentRef=cref)}, attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
    equation
      startVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, ty);
      workingArgs.maxTmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_INDEX(startVar.index+(workingArgs.numParameters+1))}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    // Modelica functions
/*
    case (DAE.CALL(path=Absyn.IDENT(ident), expLst=expList, attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, workingArgs))
    equation
      // check if function exits in workingStatargs funcNames, else append

      // process all call armugments by with expList
      // ops = collectOperationsForExpLst(ops)
      // create tmp opds and ops for contigous indices for input of modelica_call
         // operation assign_active is always created
      // create modelica_call( ident = "ident") operator
      // create operation for the modelica call with tmp opds
      (resVar, workingArgs.maxTmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, ty);
      result = OPERAND_VAR(resVar);
    then
      (inExp, (result::rest, ops, workingArgs));
*/
    case (DAE.UNARY(exp=e1), (opds, ops, workingArgs))
    equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.maxTmpIndex, Expression.typeof(e1));
      workingArgs.maxTmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_NEG(), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

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
  list<String> mathOps = {"exp", "log", "sqrt"};
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inName) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isMathFunction;

protected function isTrigMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"sin", "cos", "asin", "acos", "asinh", "acosh", "atanh"};
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inName) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isTrigMathFunction;

protected function isTransformingMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"tan","sinh","cosh","tanh"};
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inName) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isTransformingMathFunction;

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

    case OPERAND_TIME()
    then "TIME";

  end match;
end printOperandStr;

public function printOperatorStr
  input MathOperator inOp;
  output String outString;
algorithm
  outString := match inOp
  local
    Absyn.Ident ident;
    Boolean b;
    case ASSIGN_ACTIVE()
    then "assign_a";

    case ASSIGN_PARAM()
    then "assign_p";

    case ASSIGN_PASSIVE()
    then "assign_d";

    case PLUS(b)
    then "plus" + (if b then "_a" else "_d") + "_a";

    case MINUS(b)
    then "min" + (if b then "_a" else "_d") + "_a";

    case MUL(b)
    then "mult" + (if b then "_a" else "_d") + "_a";

    case DIV(b)
    then "div" + (if b then "_a" else "_d") + "_a";

    case POW()
    then "pow_op";

    case UNARY_NEG()
    then "neg_sign_a";

    case UNARY_CALL(ident)
    then ident + "_op";

  end match;
end printOperatorStr;


annotation(__OpenModelica_Interface="backend");
end MathOperation;
