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
  case simCode as SIMCODE(modelOperationData=modelOperationData,modelInfo=MODELINFO(varInfo=VARINFO(numParams=numParams))) then
    let text = (modelOperationData |> opData as OPERATIONDATA(name=name) => textFile(createAdolcTrace(opData,numParams), '<%name%>_aat.txt'))
    <<>>
  end match 
end generateAdolcAsciiTrace;

template createAdolcTrace(OperationData modelOperationData, Integer numParams)
::=
  match modelOperationData
    case operationData as OPERATIONDATA(maxTmpIndex=maxTmpIndex, independents=inds, dependents=deps) then
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
    let num_real_param = '{ op:set_numparam loc:<%numParams%> }'
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
    // define depenpendent
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

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2