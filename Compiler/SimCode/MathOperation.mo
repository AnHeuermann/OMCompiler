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
  record COND_ASSIGN
  end COND_ASSIGN;
  record COND_EQ_ASSIGN
  end COND_EQ_ASSIGN;
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
    String name;
  end OPERATIONDATA;
end OperationData;

protected uniontype WorkingStateArgs
  record WORKINGSTATEARGS
    SimCode.HashTableCrefToSimVar crefToSimVarHT;
    list<Absyn.Path> funcNames;
    Integer tmpIndex;
    Integer numVariables;
    Integer numParameters;
  end WORKINGSTATEARGS;
end WorkingStateArgs;


public function createOperationData
  input list<SimCode.SimEqSystem> inEquations;
  input SimCode.HashTableCrefToSimVar crefToSimVarHT;
  input SimCode.VarInfo varInfo;
  input String modelName;
  input DAE.FunctionTree functionTree;
  input list<SimCodeVar.SimVar> independents;
  input list<SimCodeVar.SimVar> dependents;
  output list<OperationData> outOperationData;
protected
  WorkingStateArgs workingArgs;
  Integer numVariables; 
  list<Integer> tmpLst;
  OperationData tmpOpData;
  list<OperationData> opDataFuncs;
algorithm
  try
    numVariables := 2*varInfo.numStateVars + varInfo.numAlgVars;
    workingArgs := WORKINGSTATEARGS(crefToSimVarHT, {}, numVariables, numVariables, varInfo.numParams);

    // create operation of the equations
    (tmpOpData,workingArgs) := createOperationEqns(inEquations, workingArgs);
    
    // create needed functions
    opDataFuncs := createOperationDataFuncs(workingArgs.funcNames, functionTree);

    tmpLst := list(var.index for var in independents);
    tmpOpData.independents := tmpLst;
    tmpLst := list(var.index for var in dependents);
    tmpOpData.dependents := tmpLst;

    tmpOpData.name := modelName;

    outOperationData := tmpOpData::opDataFuncs;
  else
    outOperationData := {};
  end try;
end createOperationData;


protected function createOperationDataFuncs
  input list<Absyn.Path> funcNames;
  input DAE.FunctionTree funcTree;
  output list<OperationData> outOperationData = {};
protected
  WorkingStateArgs workingArgs;
  list<Absyn.Path> funcList = funcNames;
  list<Absyn.Path> allFuncList = funcNames;
  SimCode.HashTableCrefToSimVar localHT;
  DAE.Function func;
  list<DAE.Element> inputVars, outputVars, protectedVars;
  list<DAE.Statement> bodyStmts;
  list<SimCodeVar.SimVar> inputSimVars, outputSimVars, protectedSimVars;
  Integer numVars;
  OperationData optData;
  list<Integer> tmpLst;
algorithm
  // get function for funcName
  // create OperationData for single func
  workingArgs := WORKINGSTATEARGS(HashTableCrefSimVar.emptyHashTable(), {}, 0, 0, 0);
  while not listEmpty(funcList) loop
    for funcName in funcList loop
      SOME(func) := DAE.AvlTreePathFunction.get(funcTree,funcName);

      inputVars := DAEUtil.getFunctionInputVars(func);
      outputVars :=  DAEUtil.getFunctionOutputVars(func);
      protectedVars := DAEUtil.getFunctionProtectedVars(func);
      bodyStmts := DAEUtil.getFunctionAlgorithmStmts(func);
      // create hashtable for inputs, outputs and protected
      localHT := HashTableCrefSimVar.emptyHashTable();

      inputSimVars := list( SimCodeUtil.makeTmpRealSimCodeVar(DAEUtil.varCref(v), BackendDAE.TMP_SIMVAR()) for v in inputVars);
      inputSimVars := SimCodeUtil.rewriteIndex(inputSimVars, 0);
      outputSimVars := list( SimCodeUtil.makeTmpRealSimCodeVar(DAEUtil.varCref(v), BackendDAE.TMP_SIMVAR()) for v in outputVars);
      outputSimVars := SimCodeUtil.rewriteIndex(outputSimVars, listLength(inputSimVars));
      protectedSimVars := list( SimCodeUtil.makeTmpRealSimCodeVar(DAEUtil.varCref(v), BackendDAE.TMP_SIMVAR()) for v in protectedVars);
      protectedSimVars := SimCodeUtil.rewriteIndex(protectedSimVars, listLength(inputSimVars)+listLength(outputSimVars));

      localHT := List.fold(inputSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      localHT := List.fold(outputSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      localHT := List.fold(protectedSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      numVars := listLength(inputSimVars)+listLength(outputSimVars)+listLength(protectedSimVars);

      workingArgs := WORKINGSTATEARGS(localHT, workingArgs.funcNames, numVars, numVars, 0);

      (optData, workingArgs) := createOperationsForFunction(bodyStmts, workingArgs, allFuncList);

      tmpLst := list(var.index for var in inputSimVars);
      optData.independents := tmpLst;
      tmpLst := list(var.index for var in outputSimVars);
      optData.dependents := tmpLst;

      optData.name := Absyn.pathString(funcName, "_");
      outOperationData := optData::outOperationData;
    end for;
    allFuncList := listAppend(workingArgs.funcNames, allFuncList);
    funcList := workingArgs.funcNames;
  end while;
end createOperationDataFuncs;

protected function createOperationsForFunction
  input list<DAE.Statement> funcBody;
  input WorkingStateArgs iworkingArgs;
  input list<Absyn.Path> origFuncList;
  output OperationData outOperationData;
  output WorkingStateArgs workingArgs = iworkingArgs;
protected
  list<Operation> operations;
algorithm
  (operations, workingArgs) := createOperationDataStmts(funcBody, workingArgs);
  operations := listReverse(operations);
  outOperationData := OPERATIONDATA(operations, workingArgs.tmpIndex, {}, {}, "");
end createOperationsForFunction;

protected function createOperationEqns
  input list<SimCode.SimEqSystem> inEquations;
  input WorkingStateArgs iworkingArgs;
  output OperationData outOperationData;
  output WorkingStateArgs workingArgs = iworkingArgs;
protected
  list<Operation> operations, tmpOps;
  constant Boolean debug = false;
  Integer maxTmpIndex = iworkingArgs.numVariables;
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
          workingArgs.tmpIndex = workingArgs.numVariables;
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
          workingArgs.tmpIndex = workingArgs.numVariables;
          (tmpOps, workingArgs) = createOperationDataStmts(statements, workingArgs);
          operations = listAppend(tmpOps, operations); 
        then ();

        else ();
      end matchcontinue;
      maxTmpIndex := intMax(maxTmpIndex,workingArgs.tmpIndex);
    end for;
    operations := listReverse(operations);
    outOperationData := OPERATIONDATA(operations, maxTmpIndex, {}, {},"");
  else
    fail();
  end try;
end createOperationEqns;


public function createOperationDataStmts
  input list<DAE.Statement> inStmts;
  input WorkingStateArgs iworkingArgs;
  output list<Operation> outOperations = {};
  output WorkingStateArgs workingArgs = iworkingArgs;
protected
  constant Boolean debug = false;
  Integer maxTmpIndex = iworkingArgs.numVariables;
algorithm
  try
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
        case DAE.STMT_ASSIGN(exp = rhs, exp1 = lhs) equation
          //operands = {};
          cref = Expression.expCref(lhs);
          simVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
          simVarOperand = OPERAND_VAR(simVar);
          workingArgs.tmpIndex = workingArgs.numVariables;
          (assignOperand, outOperations, workingArgs) = collectOperationsForExp(rhs, outOperations, workingArgs);
          if isTmpOperand(assignOperand) then
            op::rest = outOperations;
            op = replaceOperationResult(op, simVarOperand);
            outOperations = op::rest;
          else
            op = OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
            outOperations = op::outOperations;
          end if;
          //print(" ops: " + printOperationStr(op) + "\n");
        then ();
        //else ();
      end matchcontinue;
      maxTmpIndex := intMax(maxTmpIndex,workingArgs.tmpIndex);
    end for;
    workingArgs.tmpIndex := maxTmpIndex;
  else
    outOperations := {};
  end try;
end createOperationDataStmts;

protected function collectOperationsForFuncArgs
  input list<DAE.Exp> inExpLst;
  input list<Operand> inOpds;
  input list<Operation> inOps;
  input WorkingStateArgs iworkingArgs;
  output Operand firstArgument;
  output list<Operand> outOpds = inOpds;
  output list<Operation> outOps = inOps;
  output WorkingStateArgs workingArgs = iworkingArgs;
protected
  list<Operation> rest;
  Integer tmpIndex, argsIndex;
  SimCodeVar.SimVar simVar;
  Operand assignOperand, simVarOperand ;
  Operation op;
  Boolean first = true;
algorithm
  argsIndex := workingArgs.tmpIndex;
  workingArgs.tmpIndex := workingArgs.tmpIndex+listLength(inExpLst);

  for exp in inExpLst loop
    assignOperand::outOpds := outOpds;
    (simVar, argsIndex) := createSimTmpVar(argsIndex, Expression.typeof(exp));
    simVarOperand := OPERAND_VAR(simVar);
    op := OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
    outOps := op::outOps;
	  if first then
	    firstArgument := simVarOperand;
	    first := false;
	  end if;
  end for;
end collectOperationsForFuncArgs;

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
  //print(" traverse Exp : " + ExpressionDump.printExpStr(inExp) +"\n");
  (outExp, outTpl) := matchcontinue (inExp, inTpl)
    local
      DAE.ComponentRef cref;
      list<DAE.ComponentRef> crefList;
      SimCode.HashTableCrefToSimVar crefToSimVarHT;
      SimCodeVar.SimVar resVar, timeVar, paramVar, extraVar, startVar;
      list<Operand> opds, opdList, rest;
      list<Operation> ops;
      Operation operation;
      Operand result, firstArg;
      DAE.Exp e1,e2,e3;
      list<DAE.Exp> expList;
      DAE.Type ty;
      DAE.Operator op;
      Integer tmpIndex;
      Operand opd1, opd2, opd3;
      Absyn.Ident ident;
      String str;
      WorkingStateArgs workingArgs;
      Absyn.Path path;
      list<DAE.Var> varLst;

      constant Boolean debug = true;

    case (e1 as DAE.RCONST(), (opds, ops, workingArgs)) equation
      opds = OPERAND_CONST(e1)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (e1 as DAE.ICONST(), (opds, ops, workingArgs)) equation
      opds = OPERAND_CONST(e1)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CREF(componentRef=DAE.CREF_IDENT(ident="time"), ty=ty), (opds, ops, workingArgs)) equation
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      operation = OPERATION({OPERAND_TIME()}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CREF(componentRef=cref, ty=ty), (opds, ops, workingArgs)) equation
      paramVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      //guard
      BackendDAE.PARAM() = paramVar.varKind;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_INDEX(paramVar.index+1)}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CREF(componentRef=cref), (opds, ops, workingArgs)) equation
      resVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      opds = OPERAND_VAR(resVar)::opds;
    then
      (inExp, (opds, ops, workingArgs));

    // records
    case (DAE.CREF(componentRef=cref, ty=ty as DAE.T_COMPLEX(complexClassType=ClassInf.RECORD())), (opds, ops, workingArgs)) equation
      print("Start record case: " + ExpressionDump.printExpStr(inExp) + "\n");
      expList = Expression.expandExpression(inExp);
      crefList = list(Expression.expCref(e) for e in expList);
      opdList = list(OPERAND_VAR(BaseHashTable.get(cr, workingArgs.crefToSimVarHT)) for cr in crefList);
      print("Generated opds for record case: " + printOperandListStr(opdList) + "\n");
      opds = listAppend(opdList,opds);
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
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.RELATION(operator = op), (opds, ops, workingArgs)) equation
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createRelationOperation(op, {opd1,opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.IFEXP(expCond=e1, expThen=e2, expElse=e3), (opds, ops, workingArgs)) equation
      opd3::opd2::opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, Expression.typeof(e2));
      workingArgs.tmpIndex = tmpIndex;
      operation = OPERATION({opd1,opd2,opd3}, checkRelation(e1), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::rest;
    then
      (inExp, (opds, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT("DIVISION"), attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, workingArgs)) equation
      op = DAE.DIV(ty);
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT("atan"), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
    equation
      opd1::rest = opds;
      // tmp = operation x^2
      op = DAE.MUL(ty);
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd1}, workingArgs.tmpIndex);
      ops = operation::ops;
      // tmp = operation 1 + tmp
      op = DAE.ADD(ty);
      opd2 = OPERAND_CONST(DAE.RCONST(1.0));
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd2,result}, workingArgs.tmpIndex);
      ops = operation::ops;
      // opd2 = 1 / tmp
      op = DAE.DIV(ty);
      (operation, opd2, tmpIndex) = createBinaryOperation(op, {opd2,result}, workingArgs.tmpIndex);
      ops = operation::ops;

      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1, opd2}, UNARY_CALL("atan"), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
      guard isMathFunction(ident)
    equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
      guard isTrigMathFunction(ident)
    equation
      opd1::rest = opds;
      (extraVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1,OPERAND_VAR(extraVar)}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    case (DAE.CALL(path=Absyn.IDENT("$_start"),  expLst={DAE.CREF(componentRef=cref)}, attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
    equation
      _::rest = opds;
      startVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_INDEX(startVar.index+(workingArgs.numParameters+1))}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::rest;
    then
      (inExp, (opds, ops, workingArgs));

    // Modelica functions
    case (DAE.CALL(path=path, expLst=expList, attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, workingArgs))
    equation
      ident = Absyn.pathString(path, "_");
      // check if function exits in workingStatargs funcNames, else append
      if not List.isMemberOnTrue(path, workingArgs.funcNames, Absyn.pathEqual) then
        workingArgs.funcNames = path::workingArgs.funcNames;
      end if;

      // process all call armugments by with expList
      (firstArg, opds, ops, workingArgs) = collectOperationsForFuncArgs(expList, opds, ops, workingArgs);

      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({firstArg, OPERAND_INDEX(listLength(expList)), OPERAND_INDEX(1)}, MODELICA_CALL(ident), result);
      ops = operation::ops;
    then
      (inExp, (result::opds, ops, workingArgs));

    case (DAE.UNARY(exp=e1, operator = DAE.UMINUS()), (opds, ops, workingArgs))
    equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, Expression.typeof(e1));
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_NEG(), result);
      ops = operation::ops;
    then
      (inExp, (result::rest, ops, workingArgs));

    /* debug */
    case (_, _) guard debug
    equation
      print(" Dump not handled exp : " + ExpressionDump.printExpStr(inExp) +"\n");
      print(ExpressionDump.dumpExpStr(inExp,0) +"\n");
    then
      (inExp, inTpl);

    else
      (inExp, inTpl);
  end matchcontinue;
  //print("ready collectOperation\n");
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

protected function createRelationOperation
  input DAE.Operator operator;
  input list<Operand> inOpds;
  input Integer inIndex;
  output Operation op;
  output Operand result;
  output Integer outIndex;
protected
  Operand opd1, opd2;
  DAE.Operator newop;
algorithm
  if Expression.isLesseqOrLess(operator) then
    opd1::opd2::_ := inOpds;
  else
    opd2::opd1::_ := inOpds;
  end if;
  newop := DAE.SUB(Expression.typeofOp(operator));
  (op, result, outIndex) := createBinaryOperation(newop,{opd2,opd1},inIndex);
end createRelationOperation;

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

protected function checkRelation
  input DAE.Exp inExp;
  output MathOperator outOp;
algorithm
  outOp := match inExp

    case DAE.RELATION(operator=DAE.GREATEREQ(_))
    then COND_EQ_ASSIGN();

    case DAE.RELATION(operator=DAE.LESSEQ(_))
    then COND_EQ_ASSIGN();

    case DAE.RELATION(operator=DAE.GREATER(_))
    then COND_ASSIGN();

    case DAE.RELATION(operator=DAE.LESS(_))
    then COND_ASSIGN();

 end match;
end checkRelation;

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

protected function isArcTrigMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"asin", "acos", "atan", "asinh", "acosh", "atanh"};
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inName) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isArcTrigMathFunction;

protected function isTrigMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"sin", "cos"};
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
  input list<OperationData> operationData;
protected
  Integer i = 0;
algorithm
  for opData in operationData loop
    print("### OpeartionData : " + opData.name + "###\n");
    for op in opData.operations loop
      i := i+1;
      print(intString(i) + " " + printOperationStr(op) + "\n");
    end for;

  end for;
end dumpOperationData;

protected function printOperationStr
  input Operation inOp;
  output String outString;
algorithm
  outString := "result: " + printOperandStr(inOp.result) + ": ";
  outString := outString + printOperatorStr(inOp.operator) + " ? ";
  outString := outString + printOperandListStr(inOp.operands);
end printOperationStr;

protected function printOperandListStr
  input list<Operand> operands;
  output String outString = "";
algorithm
  for operand in operands loop
    outString := outString + printOperandStr(operand) + " ";
  end for;
end printOperandListStr;

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
    then "VAR(" + ComponentReference.printComponentRefStr(name) + " " + intString(index) + ")";

    case OPERAND_CONST(exp)
    then "CONST(" + ExpressionDump.printExpStr(exp) + ") ";

    case OPERAND_TIME()
    then "TIME";

    case OPERAND_INDEX(index)
    then "INDEX: " + intString(index);

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

    case MODELICA_CALL(ident)
    then "call fname:\"" + ident + "_aat.txt\"";

    case COND_ASSIGN()
    then "cond_assign";

    case COND_EQ_ASSIGN()
    then "cond_eq_assign";
  end match;
end printOperatorStr;


annotation(__OpenModelica_Interface="backend");
end MathOperation;
