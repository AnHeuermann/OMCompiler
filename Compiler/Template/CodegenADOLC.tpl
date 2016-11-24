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
                          odeEquations=odeEquations) then
    let()= System.tmpTickResetIndex(0,25) /* reset tmp index */
    // states are independent variables
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
    // states are independent variables
    let assign_ind = ""
    let &assign_ind += (vars.stateVars  |> var as SIMVAR(__) =>
        '{ op:assign_ind loc:<%index%> }'
    ;separator="\n")
    // derivates are dependent variables
    let assign_dep = ""
    let &assign_dep += (vars.derivativeVars  |> var as SIMVAR(__) =>
        '{ op:assign_dep loc:<%index%> }'
    ;separator="\n")
    
    let()= System.tmpTickResetIndex(0,28) /* reset ind index */
    let tickMax25 = System.tmpTickIndexReserve(28, System.tmpTickMaximum(25))
    
    
    let death_not = '{ op:death_not loc:0 loc:<%System.tmpTickMaximum(28)%> }'
    <<
    // allocation of used variables
    <%assign_zero%>
    // define independent
    <%assign_ind%>
    // operations
    <%%>
    // define depenpendent
    <%assign_dep%>
    // death_not
    <%death_not%>
    >>
  end match
end createAdolcText;

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2