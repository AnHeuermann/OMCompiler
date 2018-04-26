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
import SimCode;

// protected imports
protected
import BackendDump;
import ComponentReference;
import DAEUtil;
import Debug;
import Expression;
import ExpressionDump;
import TaskSystemDump;
import SimCodeUtil;
import System;
import Util;
import List;
import DAEDump;

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
  record EXT_DIFF_V
  end EXT_DIFF_V;
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
    Integer numParameters;
  end OPERATIONDATA;
end OperationData;

protected uniontype WorkingStateArgs
  record WORKINGSTATEARGS
    SimCode.HashTableCrefToSimVar crefToSimVarHT;
    list<Absyn.Path> funcNames;
    list<SimCode.NonlinearSystem> nlsSystems;
    Integer tmpIndex;
    Integer numVariables;
    Integer numParameters;
  end WORKINGSTATEARGS;
end WorkingStateArgs;


// Entry point for this module
public function createOperationData
  input list<SimCode.SimEqSystem> inEquations;
  input SimCode.HashTableCrefToSimVar crefToSimVarHT;
  input Integer numVariables;
  input list<SimCodeVar.SimVar> realParameters; 
  input String modelName;
  input DAE.FunctionTree functionTree;
  input list<SimCodeVar.SimVar> independents;
  input list<SimCodeVar.SimVar> dependents;
  output list<OperationData> outOperationData;
protected
  WorkingStateArgs workingArgs;
  list<Integer> tmpLst;
  OperationData tmpOpData;
  list<OperationData> opDataFuncs, opDataNLS;
  constant Boolean debug = false;
algorithm
  try
    workingArgs := WORKINGSTATEARGS(crefToSimVarHT, {}, {}, numVariables, numVariables, listLength(realParameters));

    if debug then
      print("# Equations: " + intString(listLength(inEquations)) + ".\n");
      SimCodeUtil.dumpSimEqSystemLst(inEquations, "\n");
      print("\n");
      print("createOperationData equations input: \n");
      print(Tpl.tplString3(TaskSystemDump.dumpEqs, inEquations, 0, false));
      print("-----\n");
    end if;

    // create operation of the equations
    (tmpOpData, workingArgs) := createOperationEqns(inEquations, workingArgs);

    if debug then
      print("\ncreateOperationData for NLS Systems: \n");
      Util.stringDelimitListPrintBuf(list(Absyn.pathString(str) for str in workingArgs.funcNames) , " ");
    end if;

    // create needed nls systems
    (opDataNLS, workingArgs) := createOperationDataNLS(workingArgs.nlsSystems, modelName, realParameters, workingArgs);

    if debug then
      print("\ncreateOperationData for functions: \n");
      Util.stringDelimitListPrintBuf(list(Absyn.pathString(str) for str in workingArgs.funcNames) , " ");
    end if;

    // create needed functions
    opDataFuncs := createOperationDataFuncs(workingArgs.funcNames, functionTree);
    
    opDataFuncs := listAppend(opDataNLS, opDataFuncs);

    tmpOpData := setInDepAndDepVars(independents, dependents, tmpOpData);

    tmpOpData.name := modelName;
    tmpOpData.numParameters := 1+listLength(realParameters);

    if debug then
      print("Created operations for model: " + modelName +
            " (Ind:" + intString(listLength(independents)) +
            " , Dep:" + intString(listLength(dependents)) + ")" +
            "\n");
    end if;

    outOperationData := tmpOpData::opDataFuncs;
  else
    outOperationData := {};
    Error.addMessage(Error.INTERNAL_ERROR, {"function createOperatonData failed."});
    fail();
  end try;
end createOperationData;

protected function setInDepAndDepVars
  input list<SimCodeVar.SimVar> independents;
  input list<SimCodeVar.SimVar> dependents;
  input OperationData inOpData;
  output OperationData outOpData = inOpData;
algorithm
  outOpData.independents := list(var.index for var in independents);
  outOpData.dependents := list(var.index for var in dependents);
end setInDepAndDepVars;

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

  constant Boolean debug = false;
algorithm
  // get function for funcName
  // create OperationData for single func
  workingArgs := WORKINGSTATEARGS(HashTableCrefSimVar.emptyHashTable(), {}, {}, 0, 0, 0);
  while not listEmpty(funcList) loop
    for funcName in funcList loop

      if debug then
        print("Create operation list for function : " + Absyn.pathString(funcName) + "\n");
      end if;

      SOME(func) := DAE.AvlTreePathFunction.get(funcTree, funcName);

      inputVars := DAEUtil.getFunctionInputVars(func);
      outputVars := DAEUtil.getFunctionOutputVars(func);
      protectedVars := DAEUtil.getFunctionProtectedVars(func);
      bodyStmts := DAEUtil.getFunctionAlgorithmStmts(func);
      bodyStmts := listReverse(bodyStmts);
      // create hashtable for inputs, outputs and protected
      localHT := HashTableCrefSimVar.emptyHashTable();

      inputSimVars := createSimVarsFromElements(inputVars);
      inputSimVars := SimCodeUtil.rewriteIndex(inputSimVars, 0);
      outputSimVars := createSimVarsFromElements(outputVars);
      outputSimVars := SimCodeUtil.rewriteIndex(outputSimVars, listLength(inputSimVars));
      protectedSimVars := createSimVarsFromElements(protectedVars);
      protectedSimVars := SimCodeUtil.rewriteIndex(protectedSimVars, listLength(inputSimVars)+listLength(outputSimVars));

      if debug then
        SimCodeUtil.dumpVarLst(inputSimVars, "inputs ("+intString(listLength(inputSimVars))+"): ");
        SimCodeUtil.dumpVarLst(outputSimVars, "outputs ("+intString(listLength(outputSimVars))+"): ");
        SimCodeUtil.dumpVarLst(protectedSimVars, "protected ("+intString(listLength(protectedSimVars))+"): ");
      end if;

      localHT := List.fold(inputSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      localHT := List.fold(outputSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      localHT := List.fold(protectedSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      numVars := listLength(inputSimVars) + listLength(outputSimVars) + listLength(protectedSimVars);

      workingArgs := WORKINGSTATEARGS(localHT, workingArgs.funcNames, {}, numVars, numVars, 0);

      (optData, workingArgs) := createOperationsForFunction(bodyStmts, workingArgs, allFuncList);

      tmpLst := list(var.index for var in inputSimVars);
      optData.independents := tmpLst;
      tmpLst := list(var.index for var in outputSimVars);
      optData.dependents := tmpLst;

      optData.name := Absyn.pathString(funcName, "_");
      
      optData.numParameters := 0;
      
      outOperationData := optData::outOperationData;
    end for;
    allFuncList := listAppend(workingArgs.funcNames, allFuncList);
    funcList := workingArgs.funcNames;
  end while;
end createOperationDataFuncs;

protected function createOperationDataNLS
  input list<SimCode.NonlinearSystem> nlsSysts;
  input String modelName;
  input list<SimCodeVar.SimVar> simVarParams;
  input WorkingStateArgs inWorkingStateArgs;
  output list<OperationData> outOperationData = {};
  output WorkingStateArgs outWorkingStateArgs = inWorkingStateArgs;
protected
  WorkingStateArgs workingArgs;
  SimCode.HashTableCrefToSimVar localHT;
  list<DAE.Exp> crefExps;
  list<SimCodeVar.SimVar> iterationSimVars, inputSimVars, resSimVars, innerSimVars;
  
  list<DAE.ComponentRef> resCrefs;

  Integer numVars;
  OperationData optData;
  list<Integer> tmpLst;

  constant Boolean debug = false;
algorithm
  // create OperationData for nls system
  for nlsSyst in nlsSysts loop

    // Torn system
    //   x   -> inputVars
    //   y1  -> inner vars
    //   y2  -> interation vars
    //   res -> residuals vars
    if debug then
      print("Create operation list for nlsSystem index: " + intString(nlsSyst.adolcIndex) + "\n");
    end if;

    resCrefs := list( ComponentReference.makeCrefIdent("$res" + intString(i), DAE.T_REAL_DEFAULT, {})  for i in 1:listLength(nlsSyst.crefs));
    // add res cref to residuals-simEqs
    if debug then
      print("Original Non-Linear System: \n");
      SimCodeUtil.dumpSimEqSystemLst(nlsSyst.eqs, "\n");
    end if;
    nlsSyst.eqs := makeAssignEqFromResidualEqs(nlsSyst.eqs, resCrefs);
    if debug then
      print("Non-Linear System after transformation: \n");
      SimCodeUtil.dumpSimEqSystemLst(nlsSyst.eqs, "\n");
    end if;

    // create Trace 1
    //   time, p, x -> parameters
    //   y2  -> independents
    //   res -> dependents

    localHT := HashTableCrefSimVar.emptyHashTable();

    // create iterations variables
    crefExps := list(Expression.crefToExp(cr) for cr in nlsSyst.crefs);
    crefExps := listReverse(crefExps);
    iterationSimVars := SimCodeUtil.createTempVarsforCrefs(crefExps, {});
    iterationSimVars := SimCodeUtil.rewriteIndex(iterationSimVars, 0);

    localHT := List.fold(iterationSimVars, SimCodeUtil.addSimVarToHashTable, localHT);

    // create inputs variables
    crefExps := list(Expression.crefToExp(cr) for cr in nlsSyst.inputCrefs);
    inputSimVars := SimCodeUtil.createTempVarsforCrefs(crefExps, {});
    inputSimVars := List.map1(inputSimVars, SimCodeUtil.setSimVarKind, BackendDAE.PARAM());
    inputSimVars := SimCodeUtil.rewriteIndex(inputSimVars, listLength(simVarParams));
    if debug then
      SimCodeUtil.dumpVarLst(inputSimVars, "inputs SimVars");
      SimCodeUtil.dumpVarLst(simVarParams, "parameters SimVars");
    end if;


    localHT := List.fold(simVarParams, SimCodeUtil.addSimVarToHashTable, localHT);
    localHT := List.fold(inputSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
      
    // create residual variables    
    crefExps := list(Expression.crefToExp(cr) for cr in resCrefs);
    crefExps := listReverse(crefExps);
    resSimVars := SimCodeUtil.createTempVarsforCrefs(crefExps, {});
    resSimVars := SimCodeUtil.rewriteIndex(resSimVars, listLength(iterationSimVars));
    if debug then
      SimCodeUtil.dumpVarLst(resSimVars, "resSimVars");
    end if;
    
    localHT := List.fold(resSimVars, SimCodeUtil.addSimVarToHashTable, localHT);

    //create inner variables
    crefExps := list(Expression.crefToExp(cr) for cr in nlsSyst.innerCrefs);
    innerSimVars := SimCodeUtil.createTempVarsforCrefs(crefExps, {});
    innerSimVars := SimCodeUtil.rewriteIndex(innerSimVars, listLength(iterationSimVars)+listLength(resSimVars));
    if debug then
      SimCodeUtil.dumpVarLst(innerSimVars, "inner SimVars");
    end if;


    localHT := List.fold(innerSimVars, SimCodeUtil.addSimVarToHashTable, localHT);

    numVars := listLength(iterationSimVars)+listLength(resSimVars)+listLength(innerSimVars);    
    workingArgs := WORKINGSTATEARGS(localHT, outWorkingStateArgs.funcNames, {}, numVars, numVars, 1+listLength(simVarParams)+listLength(inputSimVars));
        
    // create operation of the equations
    (optData, workingArgs) := createOperationEqns(nlsSyst.eqs, workingArgs);
    // set dep and indep  
    optData := setInDepAndDepVars(iterationSimVars, resSimVars, optData);
    optData.numParameters := 1+listLength(simVarParams)+listLength(inputSimVars);
    // set op data name
    optData.name := modelName + "_nls_" + intString(nlsSyst.adolcIndex) + "_1";
  
    outOperationData := optData::outOperationData;    

    // create Trace 3
    //   time, p, x -> parameters
    //   y2 -> independents
    //   y1 -> dependents

    // set dep and indep
    optData := setInDepAndDepVars(iterationSimVars, innerSimVars, optData);
    optData.numParameters := 1+listLength(simVarParams)+listLength(inputSimVars);
    // set op data name
    optData.name := modelName + "_nls_" + intString(nlsSyst.adolcIndex) + "_3";

    outOperationData := optData::outOperationData;

    // create Trace 2
    //   time, p, y2 -> parameters
    //   x      -> independents
    //   res,y1 -> dependents
    localHT := HashTableCrefSimVar.emptyHashTable();
    
    // x -> indep
    inputSimVars := List.map1(inputSimVars, SimCodeUtil.setSimVarKind, BackendDAE.VARIABLE());
    inputSimVars := SimCodeUtil.rewriteIndex(inputSimVars, 0);
    localHT := List.fold(inputSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
    
    // y -> params
    iterationSimVars := List.map1(iterationSimVars, SimCodeUtil.setSimVarKind, BackendDAE.PARAM());
    iterationSimVars := SimCodeUtil.rewriteIndex(iterationSimVars, listLength(simVarParams));

    localHT := List.fold(simVarParams, SimCodeUtil.addSimVarToHashTable, localHT);
    localHT := List.fold(iterationSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
    
    // res > dep
    resSimVars := SimCodeUtil.rewriteIndex(resSimVars, listLength(inputSimVars));
    localHT := List.fold(resSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
    // inner -> dep
    innerSimVars := SimCodeUtil.rewriteIndex(innerSimVars, listLength(inputSimVars)+listLength(resSimVars));
    localHT := List.fold(innerSimVars, SimCodeUtil.addSimVarToHashTable, localHT);
    
    numVars := listLength(inputSimVars)+listLength(resSimVars)+listLength(innerSimVars);
    workingArgs := WORKINGSTATEARGS(localHT, outWorkingStateArgs.funcNames, {}, numVars, numVars, 1+listLength(simVarParams)+listLength(iterationSimVars));
        
    // create operation of the equations
    (optData, workingArgs) := createOperationEqns(nlsSyst.eqs, workingArgs);
    // set dep and indep  
    optData := setInDepAndDepVars(inputSimVars, listAppend(resSimVars,innerSimVars), optData);
    optData.numParameters := 1+listLength(simVarParams)+listLength(iterationSimVars);
    // set op data name
    optData.name := modelName + "_nls_" + intString(nlsSyst.adolcIndex) + "_2";
  
    outOperationData := optData::outOperationData;   
  end for;
end createOperationDataNLS;

protected function makeAssignEqFromResidualEqs
  input list<SimCode.SimEqSystem> inEqns;
  input list<DAE.ComponentRef> inCrefs;
  output list<SimCode.SimEqSystem> outEqns = {};
protected
  list<DAE.ComponentRef> rest = inCrefs;
algorithm
  for eq in inEqns loop
    outEqns := match eq
      local
        Integer index;
        DAE.Exp exp;
        DAE.ElementSource source;
        DAE.ComponentRef cref;
      case SimCode.SES_RESIDUAL(index=index, exp=exp, source=source) equation
        cref::rest = rest;        
      then SimCode.SES_SIMPLE_ASSIGN(index, cref, exp, source, BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN)::outEqns;
      else eq::outEqns;
    end match;
  end for;
  outEqns := listReverse(outEqns);
end makeAssignEqFromResidualEqs;

protected function createSimVarsFromElements
  input list<DAE.Element> elmts; // DAE.VAR
  output list<SimCodeVar.SimVar> simVars = {};
protected
  list<SimCodeVar.SimVar> tmpVars;
  SimCodeVar.SimVar tmpVar;
  list<DAE.ComponentRef> tmpCrefs;
  list<DAE.Exp> tmpExpLst;
  DAE.ComponentRef cref;
  DAE.Type tp;
algorithm
  for elem in elmts loop
    _ := match elem
    case DAE.VAR() algorithm
    cref := DAEUtil.varCref(elem);
    // UGLY work-a-round get type of var by creating empty exp of dims
    (_, tp) := Expression.makeZeroExpression(elem.dims);
    cref := ComponentReference.crefSetType(cref, tp);
    //print("Create SimVars for Cref: " + ComponentReference.debugPrintComponentRefTypeStr(cref) + " tp: " + Types.unparseType(tp) +"\n");
    //print("inst dims: " + Util.stringDelimitList(list(intString(i) for i in  Expression.dimensionsList(elem.dims)), ",") + "\n");
    _ := match cref
      case (_) guard ComponentReference.isArrayElement(cref)
      algorithm
        tmpCrefs := ComponentReference.expandCref(cref, true);
        tmpExpLst := List.map(tmpCrefs, Expression.crefExp);
        tmpVars := SimCodeUtil.createArrayTempVar(cref, Expression.dimensionsList(elem.dims), tmpExpLst, {});
        simVars := listAppend(tmpVars, simVars);
      then ();
      else
      algorithm
        tmpVar := SimCodeUtil.makeTmpRealSimCodeVar(cref, BackendDAE.TMP_SIMVAR());
        simVars := tmpVar::simVars;
      then ();
    end match;
    then ();
    end match;
  end for;
  simVars := listReverse(simVars);
end createSimVarsFromElements;

protected function createOperationsForFunction
  input list<DAE.Statement> funcBody;
  input WorkingStateArgs inWorkingArgs;
  input list<Absyn.Path> origFuncList;
  output OperationData outOperationData;
  output WorkingStateArgs workingArgs = inWorkingArgs;
protected
  list<Operation> operations;
algorithm
  (operations, workingArgs) := createOperationDataStmts(funcBody, workingArgs);
  operations := listReverse(operations);
  outOperationData := OPERATIONDATA(operations, workingArgs.tmpIndex, {}, {}, "", 0);
end createOperationsForFunction;

protected function createOperationEqns
  input list<SimCode.SimEqSystem> inEquations;
  input WorkingStateArgs inWorkingArgs;
  output OperationData outOperationData;
  output WorkingStateArgs workingArgs = inWorkingArgs;
protected
  list<Operation> operations, tmpOps;
  Integer maxTmpIndex = inWorkingArgs.numVariables;
  list<DAE.Statement> statements;
  constant Boolean debug = false;
algorithm
  try
    operations := {};
    for eq in inEquations loop
      () := matchcontinue eq
      local
        Integer index, adolcIndex;
        DAE.Exp exp;
        list<DAE.Exp> expLst;
        DAE.ComponentRef cref;
        list<DAE.ComponentRef> crefs, inputCrefs, innerCrefs;
        list<Operation> rest;
        Operation op;
        list<Operand> operands;
        Operand assignOperand, simVarOperand, result;
        SimCodeVar.SimVar simVar;
        list<SimCodeVar.SimVar> vars;
        list<tuple<Integer, Integer, SimCode.SimEqSystem>> simJac;
        Integer i, nnz, indexB, nb,  indexA, indexX, row, col, nx, adolcIndex;
        list<Operand> intOpds = {};
        SimCode.NonlinearSystem nlSystem;

        // SIMPLE_ASSIGN
        case SimCode.SES_SIMPLE_ASSIGN(index=index, exp=exp, cref=cref) equation
          //print("collectOperationsForFuncArgs operation : " + ComponentReference.printComponentRefStr(cref) + " = " +  ExpressionDump.printExpStr(exp) +"\n");
          simVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
          simVarOperand = OPERAND_VAR(simVar);
          //print("collectOperationsForFuncArgs assign : ");
          workingArgs.tmpIndex = workingArgs.numVariables;
          ({assignOperand}, operations, workingArgs) = collectOperationsForExp(exp, operations, workingArgs);
          //print("Done with collectOperationsForExp\n");
          if isTmpOperand(assignOperand) then
            op::rest = operations;
            op = replaceOperationResult(op, simVarOperand);
            operations = op::rest;
          else
            op = OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
            operations = op::operations;
          end if;
          //print("collectOperationsForFuncArgs operation : " +  printOperationStr(op) +"\n");
        then ();

        case SimCode.SES_RESIDUAL() equation
           print("======================Hit a SES_SES_RESIDUAL\n");
        then ();

        // SES_ARRAY_CALL_ASSIGN
        case SimCode.SES_ARRAY_CALL_ASSIGN() equation
          print("======================Hit a SES_ARRAY_CALL_ASSIGN\n");
        then ();

        // SES_IFEQUATION
        case SimCode.SES_IFEQUATION() equation
          print("======================Hit a SES_IFEQUATION\n");
        then ();

        // SES_MIXED
        case SimCode.SES_MIXED() equation
          print("======================Hit a SES_MIXED\n");
        then ();

        // SES_WHEN
        case SimCode.SES_WHEN() equation
          print("======================Hit a SES_WHEN\n");
        then ();

        // SES_FOR_LOOP
        case SimCode.SES_FOR_LOOP() equation
          print("======================Hit a SES_FOR_LOOP\n");
        then ();

        // SES_INVERSE_ALGORITHM
        case SimCode.SES_INVERSE_ALGORITHM() equation
          print("======================Hit a SES_INVERSE_ALGORITHM\n");
        then ();

        // ALGORITHM, can also be (y, z) := adolcTest.f(x)
        case SimCode.SES_ALGORITHM(index=index, statements=statements) equation
          workingArgs.tmpIndex = workingArgs.numVariables;
          (tmpOps, workingArgs) = createOperationDataStmts(statements, workingArgs);
          operations = listAppend(tmpOps, operations); 
        then ();

        // SES_NONLINEAR
        case SimCode.SES_NONLINEAR(nlSystem= nlSystem as SimCode.NONLINEARSYSTEM(adolcIndex=adolcIndex, inputCrefs=inputCrefs, innerCrefs=innerCrefs, crefs=crefs)) algorithm
          //print("======================Hit a SES_NONLINEAR with index " + intString(adolcIndex) + "\n");
          
          
          indexX := workingArgs.tmpIndex;
          i := 0;
          for cr in listReverse(inputCrefs) loop
            simVar := BaseHashTable.get(cr, workingArgs.crefToSimVarHT);
            op := OPERATION({OPERAND_VAR(simVar)}, ASSIGN_ACTIVE(), OPERAND_INDEX(indexX+i));
            operations := op::operations;
            i := i + 1;
          end for;
          workingArgs.tmpIndex := workingArgs.tmpIndex+i;
          
          intOpds := listAppend({OPERAND_INDEX(adolcIndex),OPERAND_INDEX(0),OPERAND_INDEX(0),
                                 OPERAND_INDEX(1),OPERAND_INDEX(2),
                                 OPERAND_INDEX(listLength(inputCrefs)),OPERAND_INDEX(indexX),
                                 OPERAND_INDEX(listLength(innerCrefs)),OPERAND_INDEX(workingArgs.tmpIndex),
                                 OPERAND_INDEX(listLength(crefs)),OPERAND_INDEX(workingArgs.tmpIndex+listLength(innerCrefs)),
                                 OPERAND_INDEX(1)},
                                 intOpds);
          result := OPERAND_INDEX(2);
          op := OPERATION(intOpds, EXT_DIFF_V(), result);
          operations := op::operations;

          i := 0;
          for cr in listAppend(listReverse(innerCrefs), listReverse(crefs)) loop
            simVar := BaseHashTable.get(cr, workingArgs.crefToSimVarHT);
            op := OPERATION({OPERAND_INDEX(workingArgs.tmpIndex+i)}, ASSIGN_ACTIVE(), OPERAND_VAR(simVar));
            operations := op::operations;
            i := i + 1;
          end for;

          workingArgs.tmpIndex := indexX;

          workingArgs.nlsSystems := nlSystem::workingArgs.nlsSystems;
        then ();

        // SES_LINEAR
        case SimCode.SES_LINEAR(lSystem=SimCode.LINEARSYSTEM(vars=vars, beqs=expLst, simJac=simJac, jacobianMatrix=NONE(), adolcIndex=adolcIndex)) algorithm
          // create location for b
          nb := listLength(expLst);
          indexB := workingArgs.tmpIndex;
          workingArgs.tmpIndex := workingArgs.tmpIndex+nb;
          i := 0;

          for exp in expLst loop
            ({assignOperand}, operations, workingArgs) := collectOperationsForExp(exp, operations, workingArgs);
            simVarOperand := OPERAND_INDEX(indexB + i);
            if isTmpOperand(assignOperand) then
              op::rest := operations;
              op := replaceOperationResult(op, simVarOperand);
              operations := op::rest;
            elseif isConstOperand(assignOperand) then
              op := OPERATION({assignOperand}, ASSIGN_PASSIVE(), simVarOperand);
              operations := op::operations;
            else
              op := OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
              operations := op::operations;
            end if;
            i := i + 1;
          end for;


          // creat locations for A
          nnz := listLength(simJac);
          indexA := workingArgs.tmpIndex;
          workingArgs.tmpIndex := workingArgs.tmpIndex+nnz;
          i := 0;

          for tpl in simJac loop
            (row,col,SimCode.SES_RESIDUAL(exp=exp)) := tpl;
            intOpds := OPERAND_INDEX(col)::OPERAND_INDEX(row)::intOpds;
            ({assignOperand}, operations, workingArgs) := collectOperationsForExp(exp, operations, workingArgs);
            simVarOperand := OPERAND_INDEX(indexA + i);

            if isTmpOperand(assignOperand) then
              op::rest := operations;
              op := replaceOperationResult(op, simVarOperand);
              operations := op::rest;
            elseif isConstOperand(assignOperand) then
              op := OPERATION({assignOperand}, ASSIGN_PASSIVE(), simVarOperand);
              operations := op::operations;
            else
              op := OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
              operations := op::operations;
            end if;
            i := i + 1;
          end for;
          intOpds := listReverse(intOpds);

          // tmp var x
          nx := listLength(vars);
          indexX := workingArgs.tmpIndex;
          workingArgs.tmpIndex := workingArgs.tmpIndex+nx;

          // create operation for ext diff call
          i := listLength(intOpds);
          intOpds := OPERAND_INDEX(adolcIndex)::OPERAND_INDEX(i)::intOpds;
          intOpds := listAppend(intOpds, {OPERAND_INDEX(i),OPERAND_INDEX(2),OPERAND_INDEX(1),OPERAND_INDEX(nnz),
                               OPERAND_INDEX(indexA), OPERAND_INDEX(nb), OPERAND_INDEX(indexB),
                               OPERAND_INDEX(nx), OPERAND_INDEX(indexX), OPERAND_INDEX(2)});
          result := OPERAND_INDEX(1);
          op := OPERATION(intOpds, EXT_DIFF_V(), result);
          operations := op::operations;

          // assign x = tmpVars
          i := 0;
          for v in vars loop
            simVar := BaseHashTable.get(v.name, workingArgs.crefToSimVarHT);
            op := OPERATION({OPERAND_INDEX(indexX+i)}, ASSIGN_ACTIVE(), OPERAND_VAR(simVar));
            operations := op::operations;

            i := i + 1;
          end for;

        then ();

        else
          equation
            if debug then
              print("Warning not handled eqn: " + SimCodeUtil.simEqSystemString(eq) + "\n");
            else
              fail();
            end if;
          then ();

      end matchcontinue;
      maxTmpIndex := intMax(maxTmpIndex, workingArgs.tmpIndex);
    end for;
    operations := listReverse(operations);
    outOperationData := OPERATIONDATA(operations, maxTmpIndex, {}, {}, "", 0);
  else
    Error.addInternalError("createModelInfo failed", sourceInfo());
    fail();
  end try;
end createOperationEqns;

protected function createOperationDataStmts
  input list<DAE.Statement> inStmts;
  input WorkingStateArgs inWorkingArgs;
  output list<Operation> outOperations = {};
  output WorkingStateArgs workingArgs = inWorkingArgs;
protected
  constant Boolean debug = false;
  Integer maxTmpIndex = inWorkingArgs.numVariables;
algorithm
  try
    if debug then
      print("createOperation for statements input: \n");
    end if;
    for stmt in inStmts loop
      () := matchcontinue stmt
        local
          DAE.Exp rhs, lhs;
          list<DAE.Exp> expLst;
          list<Operation> rest;
          Operation op;
          Operand assignOperand, simVarOperand;
          list<Operand> assignOperands = {};

        // STMT_ASSIGN
        case DAE.STMT_ASSIGN(exp=rhs, exp1=lhs) equation
          simVarOperand = generateOperandForExp(lhs, workingArgs.crefToSimVarHT);
          workingArgs.tmpIndex = workingArgs.numVariables;
          ({assignOperand}, outOperations, workingArgs) = collectOperationsForExp(rhs, outOperations, workingArgs);
          if isTmpOperand(assignOperand) then
            op::rest = outOperations;
            op = replaceOperationResult(op, simVarOperand);
            outOperations = op::rest;
          else
            op = OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
            outOperations = op::outOperations;
          end if;
          // print(" ops: " + printOperationStr(op) + "\n");
        then ();

        // STMT_TUPLE_ASSIGN, e.g.: 4:   (y, z) := adolcTest.f(x);
        case(DAE.STMT_TUPLE_ASSIGN(expExpLst=expLst, exp=rhs)) algorithm
          if debug then
            print("======================Hit a SES_STMT_TUPLE_ASSIGN\n");
            print("Exp     | " + ExpressionDump.printExpStr(rhs) + "|\n");
            print("ExpType | " + printExpTypeStr(rhs) + " |\n");
            print("ExpLst  | " + ExpressionDump.printExpListStr(expLst) + " |\n");
          end if;

          // rhs -> list<assignOperand>
          (assignOperands, outOperations, workingArgs) := collectOperationsForExp(rhs, outOperations, workingArgs);
          if debug then
            print("Operand    | " + printOperandListStr(assignOperands) + " |\n");
            print("Operations | " + printOperationListStr(outOperations) + " |\n");
          end if;

          // lhs -> list<simVarOperand> and assign operation list<assignOperand> <--> list<simVarOperand>
          for exp in expLst loop
            simVarOperand := generateOperandForExp(exp, workingArgs.crefToSimVarHT);
            assignOperand::assignOperands := assignOperands;
            op := OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
            outOperations := op::outOperations;
          end for;
        then ();
/*
        case(DAE.STMT_TUPLE_ASSIGN(expExpLst=expLst, exp=rhs)) equation
          // rhs -> list<assignOperand>
          (assignOperand, outOperations, workingArgs) = collectOperationsForExp(rhs, outOperations, workingArgs);
          // lhs -> list<simVarOperand>
          for exp in expExpLst loop
            //
            cref = Expression.expCref(exp);
            simVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
            simVarOperand = OPERAND_VAR(simVar);

          end for;
          // assign operation list<assignOperand>  list<simVarOperand>
          op = OPERATION({assignOperand}, ASSIGN_ACTIVE(), simVarOperand);
        then()
*/
        else
        algorithm
          print("Warning not handled stmt: " + DAEDump.ppStatementStr(stmt) + "\n");
        then ();
      end matchcontinue;
      maxTmpIndex := intMax(maxTmpIndex, workingArgs.tmpIndex);
    end for;
    workingArgs.tmpIndex := maxTmpIndex;
  else
    outOperations := {};
  end try;
end createOperationDataStmts;


protected function generateOperandForExp
  "Generates an OPERAND_VAR for a given DAE.Exp and crefToSimVarHT."
  input DAE.Exp inExp;
  input SimCode.HashTableCrefToSimVar inCrefToSimVarHT;
  output Operand outSimVarOperand;
protected
  DAE.ComponentRef cref;
  SimCodeVar.SimVar simVar;
algorithm
  cref := Expression.expCref(inExp);
  simVar := BaseHashTable.get(cref, inCrefToSimVarHT);
  outSimVarOperand := OPERAND_VAR(simVar);
end generateOperandForExp;

protected function collectOperationsForFuncArgs
  input list<DAE.Exp> inExpLst;
  input list<Operand> inOpds;
  input list<Operation> inOps;
  input WorkingStateArgs inWorkingArgs;
  output Operand firstArgument;
  output list<Operand> outOpds = inOpds;
  output list<Operation> outOps = inOps;
  output WorkingStateArgs workingArgs = inWorkingArgs;
protected
  Integer tmpIndex, argsIndex;
  SimCodeVar.SimVar simVar;
  Operand assignOperand, simVarOperand ;
  Operation op;
  Boolean first = true;
  list<Operand> tmpOpds;
algorithm
  argsIndex := workingArgs.tmpIndex;
  workingArgs.tmpIndex := workingArgs.tmpIndex + listLength(inExpLst);
  tmpOpds := listReverse(List.firstN(outOpds, listLength(inExpLst)));
  outOpds := List.stripN(outOpds, listLength(inExpLst));


  for exp in inExpLst loop
    assignOperand::tmpOpds := tmpOpds;
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
  input WorkingStateArgs inWorkingArgs;
  output List<Operand> outOperands;
  output list<Operation> outOps;
  output WorkingStateArgs outWorkingArgs;
algorithm
  (_, (outOperands, outOps, outWorkingArgs)) := Expression.traverseExpBottomUp(inExp, collectOperation, ({}, inOps, inWorkingArgs));
end collectOperationsForExp;

protected function collectOperation
  input DAE.Exp inExp;
  input tuple<list<Operand>, list<Operation>, WorkingStateArgs> inTpl;
  output DAE.Exp outExp;
  output tuple<list<Operand>, list<Operation>, WorkingStateArgs> outTpl;
algorithm
  (outExp, outTpl) := matchcontinue (inExp, inTpl)
    local
      DAE.ComponentRef cref;
      list<DAE.ComponentRef> crefList;
      SimCode.HashTableCrefToSimVar crefToSimVarHT;
      SimCodeVar.SimVar resVar, timeVar, paramVar;
      list<Operand> opds, opdList, rest, results;
      list<Operation> ops, tmpOps;
      Operation operation;
      Operand result, firstArg, opd1, opd2, opd3;
      DAE.Exp e1, e2, e3;
      list<DAE.Exp> expList;
      list<DAE.Var> varLst;
      list<DAE.Subscript> subs;
      DAE.Type ty;
      DAE.Operator op;
      Integer tmpIndex;
      Absyn.Ident ident;
      Absyn.Path path;
      WorkingStateArgs workingArgs;
      constant Boolean debug = false;

    // ICONST
    case (e1 as DAE.ICONST(), (opds, ops, workingArgs)) equation
      opds = OPERAND_CONST(e1)::opds;
    then (inExp, (opds, ops, workingArgs));

    // RCONST
    case (e1 as DAE.RCONST(), (opds, ops, workingArgs)) equation
      opds = OPERAND_CONST(e1)::opds;
    then (inExp, (opds, ops, workingArgs));

    // CREF, time
    case (DAE.CREF(componentRef=DAE.CREF_IDENT(ident="time"), ty=ty), (opds, ops, workingArgs)) equation
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_TIME()}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then (inExp, (opds, ops, workingArgs));

    // CALL, start
    case (DAE.CREF(componentRef=DAE.CREF_QUAL(ident="$START",componentRef=cref),ty=ty), (opds, ops, workingArgs))
    equation
      opds = List.stripN(opds, listLength(ComponentReference.crefSubs(cref)));
      // Start var
      paramVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_INDEX(paramVar.index+(workingArgs.numParameters+1))}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then (inExp, (opds, ops, workingArgs));

    // CREF, ty
    case (DAE.CREF(componentRef=cref, ty=ty), (opds, ops, workingArgs)) equation
      paramVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      //guard
      BackendDAE.PARAM() = paramVar.varKind;
      opds = List.stripN(opds, listLength(ComponentReference.crefSubs(cref)));
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      operation = OPERATION({OPERAND_INDEX(paramVar.index+1)}, ASSIGN_PARAM(), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::opds;
    then (inExp, (opds, ops, workingArgs));

    // CREF
    case (DAE.CREF(componentRef=cref), (opds, ops, workingArgs)) equation
      resVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      opds = List.stripN(opds, listLength(ComponentReference.crefSubs(cref)));
      opds = OPERAND_VAR(resVar)::opds;
      // debug
      //print("Start cref case: " + ExpressionDump.printExpStr(inExp) + "\n");
      //print("Subs: " + ComponentReference.printComponentRef2Str("test", subs) + "\n");
      //print("res Var cref: " + printOperandListStr(opds) + "\n");
    then (inExp, (opds, ops, workingArgs));

    // CREF, records
    case (DAE.CREF(componentRef=cref), (opds, ops, workingArgs))
      guard Expression.isRecordType(ComponentReference.crefType(cref))
    equation
      //print("Start record case: " + ExpressionDump.printExpStr(inExp) + "\n");
      expList = Expression.expandExpression(inExp);
      crefList = list(Expression.expCref(e) for e in expList);
      opdList = list(OPERAND_VAR(BaseHashTable.get(cr, workingArgs.crefToSimVarHT)) for cr in crefList);
      //print("Generated opds for record case: " + printOperandListStr(opdList) + "\n");
      opds = listAppend(opdList, opds);
    then fail();
      //(inExp, (opds, ops, workingArgs));

    // CALL, der
    case (DAE.CALL(path=Absyn.IDENT("der")), (opds, ops, workingArgs)) equation
      OPERAND_VAR(resVar)::rest = opds;
      cref = ComponentReference.crefPrefixDer(resVar.name);
      resVar = BaseHashTable.get(cref, workingArgs.crefToSimVarHT);
      opds = List.stripN(opds, listLength(ComponentReference.crefSubs(cref)));
      opds = OPERAND_VAR(resVar)::rest;
    then (inExp, (opds, ops, workingArgs));

    // cref array: the list is done already
    case (DAE.ARRAY(array=expList), (opds, ops, workingArgs)) equation
      //crefList = list(Expression.expCref(e) for e in expList);
      //opdList = list(OPERAND_VAR(BaseHashTable.get(cr, workingArgs.crefToSimVarHT)) for cr in crefList);
      //opds = listAppend(opdList,opds);
    then (inExp, (opds, ops, workingArgs));

    // BINARY
    case (DAE.BINARY(operator=op), (opds, ops, workingArgs)) equation
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1, opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // BINARY, POW
    case (DAE.BINARY(operator=DAE.POW(ty=ty)), (opds, ops, workingArgs)) equation
      opd2::opd1::rest = opds;
      _ = match(opd1, opd2)
        case (OPERAND_VAR(_), OPERAND_CONST(_)) equation
          (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
          result = OPERAND_VAR(resVar);
          operation = OPERATION({opd1, opd2}, UNARY_CALL("pow"), result);
          ops = operation::ops;
        then ();
        case (OPERAND_CONST(e1), OPERAND_VAR(_)) equation
          op = DAE.MUL(ty);
          (operation, result, tmpIndex) = createBinaryOperation(op, {OPERAND_CONST(logReal(e1)), opd2}, workingArgs.tmpIndex);
          ops = operation::ops;
          operation = OPERATION({result}, UNARY_CALL("exp"), result);
          ops = operation::ops;
          then ();
        case (OPERAND_VAR(_), OPERAND_VAR(_)) equation
          (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
          result = if isTmpOperand(opd1) then opd1 else OPERAND_VAR(resVar);
          operation = OPERATION({opd1}, UNARY_CALL("log"), result);
          ops = operation::ops;
          op = DAE.MUL(ty);
          (operation, opd3, tmpIndex) = createBinaryOperation(op, {result,opd2}, workingArgs.tmpIndex);
          operation = replaceOperationResult(operation, result);
          tmpIndex = workingArgs.tmpIndex;
          ops = operation::ops;
          operation = OPERATION({result}, UNARY_CALL("exp"), result);
          ops = operation::ops;
        then ();
      end match;
      workingArgs.tmpIndex = tmpIndex;
    then (inExp, (result::rest, ops, workingArgs));

    // RELATION
    case (DAE.RELATION(operator=op), (opds, ops, workingArgs)) equation
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createRelationOperation(op, {opd1,opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // LBINARY
    case (DAE.LBINARY(operator=op), (opds, ops, workingArgs)) equation
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createLBinaryOperation(op, {opd1,opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // IFEXP
    case (DAE.IFEXP(expCond=e1, expThen=e2, expElse=e3), (opds, ops, workingArgs)) equation
      opd3::opd2::opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, Expression.typeof(e2));
      workingArgs.tmpIndex = tmpIndex;
      ({opd1,opd2,opd3}, tmpOps, tmpIndex) = createTmpLogForVals({opd1,opd2,opd3}, tmpIndex);
      ops = listAppend(tmpOps, ops);
      operation = OPERATION({opd1,opd2,opd3}, checkRelation(e1), OPERAND_VAR(resVar));
      ops = operation::ops;
      opds = OPERAND_VAR(resVar)::rest;
    then (inExp, (opds, ops, workingArgs));

    // CALL, DIVISION
    case (DAE.CALL(path=Absyn.IDENT("DIVISION"), attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, workingArgs)) equation
      op = DAE.DIV(ty);
      opd2::opd1::rest = opds;
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd2}, workingArgs.tmpIndex);
      workingArgs.tmpIndex = tmpIndex;
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // CALL, tan
    case (DAE.CALL(path=Absyn.IDENT("tan"), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs)) equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      opd2 = OPERAND_VAR(resVar);
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      opd3 = OPERAND_VAR(resVar);
      (paramVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      operation = OPERATION({opd1,OPERAND_VAR(paramVar)}, UNARY_CALL("sin"), opd2);
      ops = operation::ops;
      operation = OPERATION({opd1,OPERAND_VAR(paramVar)}, UNARY_CALL("cos"), opd3);
      ops = operation::ops;
      op = DAE.DIV(ty);
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd2,opd3}, workingArgs.tmpIndex);
      ops = operation::ops;
      workingArgs.tmpIndex = tmpIndex;
    then (inExp, (result::rest, ops, workingArgs));

    // CALL, atan
    case (DAE.CALL(path=Absyn.IDENT("atan"), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs)) equation
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

      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1, opd2}, UNARY_CALL("atan"), result);
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // CALL, asin
    case (DAE.CALL(path=Absyn.IDENT("asin"), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs)) equation
      opd1::rest = opds;
      // tmp = operation x^2
      op = DAE.MUL(ty);
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd1}, workingArgs.tmpIndex);
      ops = operation::ops;
      // tmp = operation 1 - tmp
      op = DAE.SUB(ty);
      opd2 = OPERAND_CONST(DAE.RCONST(1.0));
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd2,result}, workingArgs.tmpIndex);
      ops = operation::ops;
      // tmp = sqrt(tmp)
      operation = OPERATION({result}, UNARY_CALL("sqrt"), result);
      ops = operation::ops;
      // opd2 = 1 / tmp
      op = DAE.DIV(ty);
      (operation, opd2, tmpIndex) = createBinaryOperation(op, {opd2,result}, workingArgs.tmpIndex);
      ops = operation::ops;
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1, opd2}, UNARY_CALL("asin"), result);
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // CALL, acos
    case (DAE.CALL(path=Absyn.IDENT("acos"), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs)) equation
      opd1::rest = opds;
      // tmp = operation x^2
      op = DAE.MUL(ty);
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd1,opd1}, workingArgs.tmpIndex);
      ops = operation::ops;
      // tmp = operation 1 - tmp
      op = DAE.SUB(ty);
      opd2 = OPERAND_CONST(DAE.RCONST(1.0));
      (operation, result, tmpIndex) = createBinaryOperation(op, {opd2,result}, workingArgs.tmpIndex);
      ops = operation::ops;
      // tmp = sqrt(tmp)
      operation = OPERATION({result}, UNARY_CALL("sqrt"), result);
      ops = operation::ops;
      // opd2 = - 1 / tmp
      opd2 = OPERAND_CONST(DAE.RCONST(-1.0));
      op = DAE.DIV(ty);
      (operation, opd2, tmpIndex) = createBinaryOperation(op, {opd2,result}, workingArgs.tmpIndex);
      ops = operation::ops;
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1, opd2}, UNARY_CALL("acos"), result);
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // CALL, mathFunction
    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
      guard isMathFunction(ident)
    equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // CALL, trigMathFunction
    case (DAE.CALL(path=Absyn.IDENT(ident), attr=DAE.CALL_ATTR(builtin=true, ty=ty)), (opds, ops, workingArgs))
      guard isTrigMathFunction(ident)
    equation
      opd1::rest = opds;
      (paramVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, ty);
      (resVar, tmpIndex) = createSimTmpVar(tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1,OPERAND_VAR(paramVar)}, UNARY_CALL(ident), result);
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // Modelica functions
    // CALL
    case (DAE.CALL(path=path, expLst=expList, attr=DAE.CALL_ATTR(ty=ty)), (opds, ops, workingArgs)) equation
      ident = Absyn.pathString(path, "_");
      // check if function exits in workingStatargs funcNames, else append
      if not List.isMemberOnTrue(path, workingArgs.funcNames, Absyn.pathEqual) then
        workingArgs.funcNames = path::workingArgs.funcNames;
      end if;
      //print("collectOperationsForFuncArgs pre opds : " +  printOperandListStr(opds) +"\n");

      // process all call armugments by with expList
      //print("collectOperationsForFuncArgs for : " + ExpressionDump.printExpListStr(expList) +"\n");
      (firstArg, opds, ops, workingArgs) = collectOperationsForFuncArgs(expList, opds, ops, workingArgs);
      //print("collectOperationsForFuncArgs opds : " +  printOperandListStr(opds) +"\n");

      (results, tmpIndex) = createOperandVarLst(workingArgs.tmpIndex, ty);
      workingArgs.tmpIndex = tmpIndex;
      result = List.first(results);
      opds = listAppend(results, opds);
      //print("collectOperationsForFuncArgs opds : " +  printOperandStr(result) +"\n");

      // Get number of result vars
      tmpIndex = numResults(ty);
      if debug then
        print("numRes = " + intString(tmpIndex) + ".\n");
        print("Type   = " + Types.printTypeStr(ty) + ".\n");
      end if;
      operation = OPERATION({firstArg, OPERAND_INDEX(listLength(expList)), OPERAND_INDEX(tmpIndex)}, MODELICA_CALL(ident), result);
      //print("collectOperationsForFuncArgs operation : " +  printOperationStr(operation) +"\n");
      ops = operation::ops;
    then (inExp, (opds, ops, workingArgs));

    // UNARY, minus
    case (DAE.UNARY(exp=e1, operator=DAE.UMINUS()), (opds, ops, workingArgs)) equation
      opd1::rest = opds;
      (resVar, tmpIndex) = createSimTmpVar(workingArgs.tmpIndex, Expression.typeof(e1));
      workingArgs.tmpIndex = tmpIndex;
      result = OPERAND_VAR(resVar);
      operation = OPERATION({opd1}, UNARY_NEG(), result);
      ops = operation::ops;
    then (inExp, (result::rest, ops, workingArgs));

    // debug
    case (_, _) guard debug
    equation
      print("Dump not handled exp : " + ExpressionDump.printExpStr(inExp) + "\n");
      print(ExpressionDump.dumpExpStr(inExp, 0) + "\n");
    then (inExp, inTpl);

    else
      (inExp, inTpl);
  end matchcontinue;
  //print("ready collectOperation\n");
end collectOperation;

protected function numResults
  "Todo: Is there an existing function in Frontend oder Backend?
   Todo: Implement the missing cases.
  Returns the number of result variables for DAE.CALL"
  input DAE.Type inType;
  output Integer numRes = 0;
algorithm
  numRes := match inType
    local
      list<DAE.Type> types;

    // TODO:
    // case T_REAL etc.

    case DAE.T_TUPLE(types=types) then listLength(types);
    else 1;
  end match;
end numResults;

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

    // ADD
    case DAE.ADD(ty) equation
      (resVar, outIndex) = createSimTmpVar(inIndex, ty);
      result = OPERAND_VAR(resVar);
      (opds, extraOps, isActive, isCommuted) = checkOperand(inOpds);
      op = OPERATION(opds, PLUS(isActive), result);
    then ();

    // SUB
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

    // MUL
    case DAE.MUL(ty) equation
      (resVar, outIndex) = createSimTmpVar(inIndex, ty);
      result = OPERAND_VAR(resVar);
      (opds, extraOps, isActive, isCommuted) = checkOperand(inOpds);
      op = OPERATION(opds, MUL(isActive), result);
    then ();

    // DIV
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
  (op, result, outIndex) := createBinaryOperation(newop, {opd2,opd1}, inIndex);
end createRelationOperation;

protected function createLBinaryOperation
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
  opd1::opd2::_ := inOpds;
  if Expression.operatorEqual(operator,DAE.AND(Expression.typeofOp(operator))) then
    newop := DAE.MUL(Expression.typeofOp(operator));
  elseif Expression.operatorEqual(operator,DAE.OR(Expression.typeofOp(operator))) then
    newop := DAE.ADD(Expression.typeofOp(operator));
  else
    Error.addMessage(Error.INTERNAL_ERROR,{"MathOperation.createLBinaryOperation failed for op : " + ExpressionDump.binopSymbol(operator)});
  end if;

  (op, result, outIndex) := createBinaryOperation(newop, {opd2,opd1}, inIndex);
end createLBinaryOperation;

protected function createSimTmpVar
 input Integer inIndex;
 input DAE.Type inType;
 output SimCodeVar.SimVar simvar;
 output Integer nextIndex;
protected
  DAE.ComponentRef cref;
algorithm
  cref := ComponentReference.makeCrefIdent("$tmpOpVar"+"_"+intString(inIndex), inType, {});
  simvar := SimCodeUtil.makeTmpRealSimCodeVar(cref, BackendDAE.TMP_SIMVAR(), inIndex);
  nextIndex := inIndex + 1;
end createSimTmpVar;

protected function createOperandVarLst
  "Creates a list of operands for the given DAE.type"
  input Integer inIndex;
  input DAE.Type inType;
  output list<Operand> operands = {};
  output Integer nextIndex = inIndex;
protected
  DAE.ComponentRef cref;
  Integer numVars;
  SimCodeVar.SimVar simVar;
algorithm
  numVars := numResults(inType);
  for i in 1:numVars loop
    cref := ComponentReference.makeCrefIdent("$tmpOpVar"+"_"+intString(nextIndex), inType, {});
    simVar := SimCodeUtil.makeTmpRealSimCodeVar(cref, BackendDAE.TMP_SIMVAR(), nextIndex);
    operands := OPERAND_VAR(simVar)::operands;
    nextIndex := nextIndex + 1;
  end for;
  operands := listReverse(operands);
end createOperandVarLst;

protected function createTmpLogForVals
  "Creates a list of operands for the given DAE.type"
  input list<Operand> inOpds;
  input Integer inIndex;
  output list<Operand> operands = {};
  output list<Operation> extraOps = {};
  output Integer nextIndex = inIndex;
protected
  SimCodeVar.SimVar simVar;
  DAE.Exp e;
  DAE.ComponentRef cref;
algorithm
  for opd in inOpds loop
    _ := match(opd)
      case (OPERAND_CONST(const=e)) equation
        cref = ComponentReference.makeCrefIdent("$tmpOpVar"+"_"+intString(nextIndex), Expression.typeof(e), {});
	      simVar = SimCodeUtil.makeTmpRealSimCodeVar(cref, BackendDAE.TMP_SIMVAR(), nextIndex);
	      extraOps = OPERATION({opd}, ASSIGN_PASSIVE(), OPERAND_VAR(simVar))::extraOps;
	      operands = OPERAND_VAR(simVar)::operands;
	      nextIndex = nextIndex + 1;
	    then ();
      else equation
        operands = opd::operands;
      then ();
    end match;
  end for;
  operands := listReverse(operands);
end createTmpLogForVals;

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
      outOpds = {opd2, opd1};
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

    case DAE.LBINARY(operator=DAE.AND(_))
    then COND_ASSIGN();

    case DAE.LBINARY(operator=DAE.OR(_))
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

protected function isConstOperand
  input Operand inOperand;
  output Boolean outBoolean;
algorithm
  outBoolean := matchcontinue inOperand
  local

    case OPERAND_CONST()
    then true;

    else false;
  end matchcontinue;
end isConstOperand;

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
  outBool := isKindOf(inName, mathOps);
end isMathFunction;

protected function isArcTrigMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"asin", "acos", "atan", "asinh", "acosh", "atanh"};
algorithm
  outBool := isKindOf(inName, mathOps);
end isArcTrigMathFunction;

protected function isTrigMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"sin", "cos"};
algorithm
  outBool := isKindOf(inName, mathOps);
end isTrigMathFunction;

protected function isKindOf
  "Returns true if the given operation is an element of the given set of math operations/strings."
  input String inOp;
  input list<String> mathOps;
  output Boolean outBool;
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inOp) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isKindOf;

protected function isTransformingMathFunction
  input String inName;
  output Boolean outBool;
protected
  list<String> mathOps = {"tan", "sinh", "cosh", "tanh"};
algorithm
  outBool := false;
  for op in mathOps loop
    if stringCompare(op, inName) == 0 then
      outBool := true;
      break;
    end if;
  end for;
end isTransformingMathFunction;

protected function logReal
  input DAE.Exp inExp;
  output DAE.Exp outExp;
algorithm
  outExp := match inExp
    local
      Integer i;
      Real r;

    case DAE.ICONST(i) equation
      r = log(i);
    then DAE.RCONST(r);

    case DAE.RCONST(r) equation
      r = log(r);
    then DAE.RCONST(r);

    else fail();
  end match;
end logReal;


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

protected function printOperationListStr
  input list<Operation> inOperations;
  output String outString = "";
algorithm
  outString := "OperationList----- \n";
  for operation in inOperations loop
    outString := outString + "\t" + printOperationStr(operation) + "\n";
  end for;
  outString := outString + "OperationList----- \n";
end printOperationListStr;

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

    case EXT_DIFF_V()
    then "ext_diff_v";
  end match;
end printOperatorStr;

protected function printExpTypeStr
  "MS: Copied from ExpressionDump.mo (protected!) for debugging.
   Prints out the name of the expression uniontype to a string."
  input DAE.Exp inExp;
  output String outString;
algorithm
  outString := match(inExp)
    case DAE.ICONST(_) then "ICONST";
    case DAE.RCONST(_) then "RCONST";
    case DAE.SCONST(_) then "SCONST";
    case DAE.BCONST(_) then "BCONST";
    case DAE.ENUM_LITERAL() then "ENUM_LITERAL";
    case DAE.CREF() then "CREF";
    case DAE.BINARY() then "BINARY";
    case DAE.UNARY() then "UNARY";
    case DAE.LBINARY() then "LBINARY";
    case DAE.LUNARY() then "LUNARY";
    case DAE.RELATION() then "RELATION";
    case DAE.IFEXP() then "IFEXP";
    case DAE.CALL() then "CALL";
    case DAE.PARTEVALFUNCTION() then "PARTEVALFUNCTION";
    case DAE.ARRAY() then "ARRAY";
    case DAE.MATRIX() then "MATRIX";
    case DAE.RANGE() then "RANGE";
    case DAE.TUPLE() then "TUPLE";
    case DAE.CAST() then "CAST";
    case DAE.ASUB() then "ASUB";
    case DAE.TSUB() then "TSUB";
    case DAE.SIZE() then "SIZE";
    case DAE.CODE() then "CODE";
    case DAE.EMPTY() then "EMPTY";
    case DAE.REDUCTION() then "REDUCTION";
    case DAE.LIST() then "LIST";
    case DAE.CONS() then "CAR";
    case DAE.META_TUPLE() then "META_TUPLE";
    case DAE.META_OPTION() then "META_OPTION";
    case DAE.METARECORDCALL() then "METARECORDCALL";
    case DAE.MATCHEXPRESSION() then "MATCHEXPRESSION";
    case DAE.BOX() then "BOX";
    case DAE.UNBOX() then "UNBOX";
    case DAE.SHARED_LITERAL() then "SHARED_LITERAL";
    case DAE.PATTERN() then "PATTERN";
    else "#UNKNOWN EXPRESSION#";
  end match;
end printExpTypeStr;


annotation(__OpenModelica_Interface="backend");
end MathOperation;
