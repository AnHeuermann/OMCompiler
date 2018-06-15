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
import ExpressionDumpTpl;

/* public */ template generateAdolcAsciiTrace(SimCode simCode)
  "Generates ADOL-C ascii trace file"
::=
  match simCode
  case simCode as SIMCODE(modelOperationData=modelOperationData) then
    let text = (modelOperationData |> opData as OPERATIONDATA(name=name) => textFile(createAdolcTrace(opData), '<%name%>_aat.txt'))
    let text2 = (modelOperationData |> opData as OPERATIONDATA(__) => createPatternFiles(opData))
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

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2