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
  case simCode as SIMCODE(__) then
    let text = createAdolcText(simCode)
    let()=textFile(text, '<%fileNamePrefix%>_adolcAsciiTrace.txt')
    ""
  end match 
end generateAdolcAsciiTrace;

template createAdolcText(SimCode simCode)
::=
  match simCode
  case simCode as SIMCODE(modelInfo=MODELINFO(vars=vars as SIMVARS(__),
                                              varInfo=varInfo as VARINFO(__)),
                          modelOperationData=modelOperationData) then
    let maxTmpIndex = match modelOperationData
                      case SOME(operationData as
                                OPERATIONDATA(maxTmpIndex=maxTmpIndex, independents=inds, dependents=deps)) then
                         '<%maxTmpIndex%>'
                      end match
    //let()= System.tmpTickResetIndex(0,25) /* reset tmp index */
    // states are independent variables
    /*
    let assign_zero = ""
    let &assign_zero += (vars.stateVars |> var as  SIMVAR(__) =>
            '{ op:assign_d_zero loc:<%index%> }'
    ;separator="\n")
    let &assign_zero += "\n"
    let &assign_zero += (vars.derivativeVars  |> var as SIMVAR(__) =>
        '{ op:assign_d_zero loc:<%index%> }'
    ;separator="\n")
    let &assign_zero += "\n"
    let &assign_zero += (vars.algVars  |> var as SIMVAR(__) =>
        '{ op:assign_d_zero loc:<%index%> }'
    ;separator="\n")
    */
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
    let num_real_param = '{ op:set_numparam loc:<%varInfo.numParams%> }'
    let operations = match modelOperationData
                      case SOME(operationData as
                                OPERATIONDATA(operations=operations)) then
                      let opsText = ""
                      let &opsText += (operations |> op as OPERATION(__)
                                       => createOperatorText(op)
                                       ;separator="\n")
                      '<%opsText%>'
                    end match
    <<
    // allocation of used variables
    <%assign_zero%>
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
end createAdolcText;

template createOperatorText(MathOperation.Operation op)
::= match op
    case OPERATION(operator=operator,operands=operands,result=result) then
       match operator
         case ASSIGN_PARAM() then
             let operStr = 'assign_p'
             let locsStr = ""
             let &locsStr += (operands |> opd
                      as OPERAND_VAR(variable=variable as SimCodeVar.SIMVAR(index=index)) => ' loc:<%intAdd(index,1)%>')
             let &locsStr += (operands |> opd as OPERAND_TIME() => ' loc:0')
             let &locsStr += match result
                        case OPERAND_VAR(variable=variable as
                                     SimCodeVar.SIMVAR(index=index)) then
                        ' loc:<%index%>'
                        end match
             <<{ op:<%operStr%><%locsStr%> }>>

         else
    let operStr = MathOperation.printOperatorStr(operator)
    let locsStr = ""
    let &locsStr += (operands |> opd as OPERAND_VAR(variable=variable as SimCodeVar.SIMVAR(index=index)) => 'loc:<%index%> ')
    let &locsStr += match result
                    case OPERAND_VAR(variable=variable as
                                     SimCodeVar.SIMVAR(index=index)) then
                    'loc:<%index%>'
                    end match
    let valStr = ""
    let &valStr += (operands |> opd as OPERAND_CONST(const=const)
                    => 'val:<%ExpressionDumpTpl.dumpExp(const,"")%> ')
    <<{ op:<%operStr%> <%locsStr%> <%valStr%>}>>
    else
    <<>>
    end match
end createOperatorText;

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2