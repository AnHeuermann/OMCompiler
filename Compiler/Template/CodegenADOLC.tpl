// This file defines templates for transforming Modelica/MetaModelica code to C
// code. They are used in the code generator phase of the compiler to write
// target code.
//
// CodegenC.tpl has the root template translateModel while
// this template contains only translateFunctions.
// These templates do not return any
// result but instead write the result to files. All other templates return
// text and are used by the root templates (most of them indirectly).

package CodegenADOLC

import interface SimCodeTV;
import CodegenUtil.*;
import CodegenCFunctions.*;
import ExpressionDumpTpl;

/* public */ template generateAdolcAsciiTrace(SimCode simCode)
  "Generates ADOL-C ascii trace file"
::=
  match simCode
  case simCode as SIMCODE(modelOperationData=modelOperationData, fileNamePrefix=fileNamePrefix) then
    let _ = (modelOperationData |> opData as OPERATIONDATA(name=name) => textFile(createAdolcTrace(opData), '<%name%>_aat.txt'))
    let _ = (modelOperationData |> opData as OPERATIONDATA(__) => createPatternFiles(opData))
    let _ =  match List.first(modelOperationData)
                case OPERATIONDATA(extFuncNames=extFuncNames) then
                    createExtFunctionsFile(fileNamePrefix, extFuncNames)
    <<>>
  end match 
end generateAdolcAsciiTrace;

template createAdolcTrace(OperationData modelOperationData)
::=
  match modelOperationData
    case operationData as OPERATIONDATA(maxTmpIndex=maxTmpIndex, independents=inds, dependents=deps, numRealParameters=numRealParameters) then
    let tmpIndex='<%maxTmpIndex%>'
    // states are independent variables
    let assign_ind = ""
    let &assign_ind += (inds  |> i as Integer =>
        '{ op:assign_ind loc:<%i%> }'
    ;separator="\n")
    // derivates are dependent variables
    let assign_dep = ""
    let &assign_dep += (deps  |> i as Integer =>
        '{ op:assign_dep loc:<%i%> }'
    ;separator="\n")

    let death_not = '{ op:death_not loc:0 loc:<%maxTmpIndex%> }'
    let num_real_param = '{ op:set_numparam loc:<%numRealParameters%> }'
    let operations = match modelOperationData
                      case operationData as
                                OPERATIONDATA(operations=operations) then
                      let opsText = ""
                      let &opsText += (operations |> op as OPERATION(__)
                                       => createOperatorText(op)
                                       ;separator="\n")
                      '<%opsText%>'
                    end match
    <<
    // allocation of used variables
    // define independent
    <%assign_ind%>
    // operations
    <%operations%>
    // define dependent
    <%assign_dep%>
    // death_not
    <%death_not%>
    // num real parameters
    <%num_real_param%>
    >>
  end match
end createAdolcTrace;

template createOperatorText(MathOperation.Operation op)
::= match op
    case OPERATION(operator=operator,operands=operands,result=result) then
    let operStr = MathOperation.printOperatorStr(operator)
    let locsStr = ""
    let &locsStr += (operands |> opd => createOperandText(opd))
    let &locsStr += createOperandText(result)
    let valStr = ""
    let &valStr += (operands |> opd as OPERAND_CONST(const=const)
                    => 'val:<%ExpressionDumpTpl.dumpExp(const,"")%> ')
    <<{ op:<%operStr%> <%locsStr%> <%valStr%>}>>
    else
    <<>>
    end match
end createOperatorText;

template createOperandText(Operand opd)
::=match opd
     case OPERAND_VAR(variable=variable as SimCodeVar.SIMVAR(index=index)) then
       << loc:<%index%>>>
     case OPERAND_INDEX(i=index) then
       << loc:<%index%>>>
     case OPERAND_TIME() then
       << loc:0>>
     else
       <<>>
  end match

end createOperandText;

template createPatternFiles(OperationData modelOperationData)
::=match modelOperationData
    case OPERATIONDATA(name=name,linSysPat=linSysPattern) then
    let text = (linSysPattern |> lsPat as LINSYSPATTERN(adolcIndex=adolcIndex)
                => textFile(createOnePattern(lsPat), '<%name%>_ls_<%adolcIndex%>_pat.txt'))
    <<>>
  end match
end createPatternFiles;

template createOnePattern(LinSysPattern linSysPat)
::=match linSysPat
    case LINSYSPATTERN(adolcIndex=adolcIndex,pattern=pattern) then
    let text = (pattern |> n as Integer => '<%n%>';separator=" ")
    << <%text%> >>
  end match
end createOnePattern;

template createExtFunctionsFile(String fileNamePrefix, list<tuple<SimCodeFunction.Function, tuple<list<ArgsIndices>, SimCodeFunction.Function>, Integer>> extFuncNames)
::=
match extFuncNames
case {} then
textFile("extern \"C\" void init_modelica_external_functions() {}", '<%fileNamePrefix%>_extFuncs.cpp')
else
textFile(createExternalFunctionCalls(fileNamePrefix, extFuncNames), '<%fileNamePrefix%>_extFuncs.cpp')
end createExtFunctionsFile;

template createExternalFunctionCalls(String fileNamePrefix, list<tuple<SimCodeFunction.Function, tuple<list<ArgsIndices>, SimCodeFunction.Function>, Integer>> extFuncNames)
::=

let funcClassDecl = ( extFuncNames |> funcNameTuple => createExtFuncClass(funcNameTuple); separator="\n\n")
let extFuncInit = (extFuncNames |> funcNameTuple => createExtFuncInit(funcNameTuple); separator="\n")
<<
/* External function calls for adolc */
#include "util/modelica.h"

#include "<%fileNamePrefix%>_includes.h"
#include "<%fileNamePrefix%>_functions.h"


#include <adolc/edfclasses.h>

#include <forward_list>

using std::forward_list;

class ModelicaExtFunc : public EDFobject {
public:
    ModelicaExtFunc() : EDFobject(){  }
    virtual ~ModelicaExtFunc() {}
    virtual int function(int n, double *x, int m, double *y) { return this->zos_forward(n,x,m,y); }
    virtual int zos_forward(int n, double *x, int m, double *y) = 0;
    virtual int fos_forward(int n, double *dp_x, double *dp_X, int m, double *dp_y, double *dp_Y)  = 0;
    virtual int fov_forward(int n, double *dp_x, int p, double **dpp_X, int m, double *dp_y, double **dpp_Y);
    virtual int fos_reverse(int m, double *dp_U, int n, double *dp_Z, double *dp_x, double *dp_y) {return 0;}
    virtual int fov_reverse(int m, int p, double **dpp_U, int n, double **dpp_Z, double *dp_x, double *dp_y) {return 0;}
};

int ModelicaExtFunc::fov_forward(int n, double *dp_x, int p, double **dpp_X, int m, double *dp_y, double **dpp_Y) {
    int ret;
    double *dp_X, *dp_Y;
    dp_X = (double*)calloc(n,sizeof(double));
    dp_Y = (double*)calloc(m,sizeof(double));
    for(int i =0; i<p; i++) {
         for (int j=0; j<n; j++)
             dp_X[j]=dpp_X[j][i];
         int iret = fos_forward(n,dp_x,dpp_X[i],m,dp_y,dp_Y);
         for (int j=0; j<m; j++)
             dpp_Y[j][i] = dp_Y[j];
         MINDEC(ret,iret);
    }
    return ret;
}

static forward_list<ModelicaExtFunc*> extfunclist;

<%funcClassDecl%>

extern "C" void init_modelica_external_functions() {
    ModelicaExtFunc* fp;
    <%extFuncInit%>
}
>>
end createExternalFunctionCalls;

template createExtFuncInit(tuple<SimCodeFunction.Function, tuple<list<ArgsIndices>, SimCodeFunction.Function>, Integer> funcTuple)
::=
match funcTuple
case (f as SimCodeFunction.EXTERNAL_FUNCTION(name=name,extArgs=extArgs, extReturn=extReturn),
      (derInputExtArgs,
        derf as  SimCodeFunction.EXTERNAL_FUNCTION(name=derName, extArgs= derExtArgs, extReturn=derExtReturn))
      , index
      ) then
let nameStr = underscorePath(name)
<<
fp = new MEF_<%nameStr%>;
extfunclist.push_front(fp);
assert(fp->get_index() == <%index%>);
>>
end createExtFuncInit;

template createExtFuncClass(tuple<SimCodeFunction.Function, tuple<list<ArgsIndices>, SimCodeFunction.Function>, Integer> funcTuple)
::=
match funcTuple
case (f as SimCodeFunction.EXTERNAL_FUNCTION(name=name,extArgs=extArgs, extReturn=extReturn),
      (derInputExtArgs,
        derf as  SimCodeFunction.EXTERNAL_FUNCTION(name=derName, extArgs= derExtArgs, extReturn=derExtReturn))
      , _
      ) then
      
let inputs = (extArgs |> arg hasindex i0 =>  match arg 
                    case SimCodeFunction.SIMEXTARG(type_=T_COMPLEX(complexClassType=EXTERNAL_OBJ(__))) then
                    if isInput then
                        let typeStr = extType(type_, isInput, isArray)
                        'memcpy(&<%extVarName(cref)%>,&dp_x[<%i0%>],sizeof(<%typeStr%>));'
                    else
                        ''
                    case SimCodeFunction.SIMEXTARG(__) then
                    if isInput then 
                        let typeStr = extType(type_, isInput, isArray)
                        '<%extVarName(cref)%> = (<%typeStr%>) dp_x[<%i0%>];'
                    else 
                        ''
                    ; separator="\n")

let output = match extReturn
               case SimCodeFunction.SIMEXTARG(__) then
                 'dp_y[0] = (<%extReturnType(extReturn)%>) <%extVarName(cref)%>;'

let extraOuts = (extArgs |> arg =>  match arg 
                    case sExt as SimCodeFunction.SIMEXTARG(__) then
                    match outputIndex
                    case 0 then ''
                    case _ then 
                    let typeStr = extType(sExt.type_, sExt.isInput, sExt.isArray)
                    'dp_y[<%outputIndex%>] = (<%typeStr%>) <%extVarName(sExt.cref)%>;'
                    ; separator="\n") 
let preExp = ""
let varDecls = ""
let varInit = ""
let auxFuncs = "" 

let extFuncCall = extFunCall(f, &preExp, &varDecls, &varInit, &auxFuncs)
let returnDecl = match extReturn
                    case SimCodeFunction.SIMEXTARG(__) then
                        '<%expTypeFlag(type_, 2)%> _<%crefStr(cref)%>;'

let nameStr = underscorePath(name)

let derInputs = (derInputExtArgs |> derArg  =>
                  match derArg
                  case arg as MathOperation.ARGS_INDICES(__) then
                    match arg.argument
                    case SimCodeFunction.SIMEXTARG(type_=T_COMPLEX(complexClassType=EXTERNAL_OBJ(__))) then
                    if isInput then
                        let typeStr = extType(type_, isInput, isArray)
                        '<%extVarName(cref)%> = reinterpret_cast< <%typeStr%> > (static_cast<size_t> (dp_X[<%arg.index%>]));'
                    case sExt as SimCodeFunction.SIMEXTARG(__) then
                    if isInput then
                      let typeStr = extType(sExt.type_, sExt.isInput, sExt.isArray)
                      '<%extVarName(sExt.cref)%> = (<%typeStr%>) dp_X[<%arg.index%>];'
                ; separator="\n")

let derOutput = match derExtReturn
               case SimCodeFunction.SIMEXTARG(__) then
                 'dp_Y[0] = (<%extReturnType(derExtReturn)%>) <%extVarName(cref)%>;'

let derExtraOuts = (derExtArgs |> arg =>  match arg 
                    case sExt as SimCodeFunction.SIMEXTARG(__) then
                    match outputIndex
                    case 0 then ''
                    case _ then 
                    let typeStr = extType(sExt.type_, sExt.isInput, sExt.isArray)
                    'dp_Y[<%outputIndex%>] = (<%typeStr%>) <%extVarName(sExt.cref)%>;'
                    ; separator="\n") 
let derPreExp = ""
let derVarDecls = ""
let derVarInit = ""
let derAuxFuncs = "" 

let derExtFuncCall = extFunCall(derf, &derPreExp, &derVarDecls, &derVarInit, &derAuxFuncs)
let derReturnDecl = match derExtReturn
                    case SimCodeFunction.SIMEXTARG(__) then
                        '<%expTypeFlag(type_, 2)%> _<%crefStr(cref)%>;'

<<
class MEF_<%nameStr%>: public ModelicaExtFunc {
public :
    MEF_<%nameStr%>() : ModelicaExtFunc() {}
    virtual ~MEF_<%nameStr%>() {}
    virtual int zos_forward(int n, double *dp_x, int m, double *dp_y);
    virtual int fos_forward(int n, double *dp_x, double *dp_X, int m, double *dp_y, double *dp_Y);
};

int MEF_<%nameStr%>::zos_forward(int n, double *dp_x, int m, double *dp_y) {

/* varDecls */
<%varDecls%>
/* varInit */
<%varInit%>
/* auxFuncs */
<%auxFuncs%>
/* return decl */
<%returnDecl%>

/* inputs => preExp */
<%inputs%>


/* external function call */
<%extFuncCall%>

/* get outputs */
<%output%>
<%extraOuts%>

return 0;
}

int MEF_<%nameStr%>::fos_forward(int n, double *dp_x, double *dp_X, int m, double *dp_y, double *dp_Y) {
{ /* extra block for local variables */
/* varDecls */
<%varDecls%>
/* varInit */
<%varInit%>
/* auxFuncs */
<%auxFuncs%>
/* return decl */
<%returnDecl%>

/* inputs => preExp */
<%inputs%>
/* external function call */
<%extFuncCall%>
/* get outputs */
<%output%>
<%extraOuts%>
} /* extra block end */

{ /* extra block for local variables */
/* derVarDecls */
<%derVarDecls%>
/* dervarInit */
<%derVarInit%>
/* derAuxFuncs */
<%derAuxFuncs%>
/* der return decl */
<%derReturnDecl%>
/* inputs => preExp */
<%inputs%>
/* derInputs */
<%derInputs%>

/* der external function call */
<%derExtFuncCall%>

/* get derOutputs */
<%derOutput%>
<%derExtraOuts%>
} /* extra block end */
return 0;
}
>>
end createExtFuncClass;

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2