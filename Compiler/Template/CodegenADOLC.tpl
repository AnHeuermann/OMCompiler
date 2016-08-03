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
  case simCode as SIMCODE(varInfo=VARINFO(numStateVars=numStateVars, numAlgVars=numAlgVars),
                          modelInfo=MODELINFO(vars=SIMVARS(__)),
                          odeEquations=odeEquations) then
    let()= System.tmpTickResetIndex(0,25) /* reset tmp index */
    // states are independent variables
    let assign_zero = (vars.stateVars  |> var =>
        match var
          case SIMVAR(__) then
            '{ op:assign_d_zero loc:<%index%> }'
    ;separator="\n")
    let assign_zero += (vars.derivativeVars  |> var =>
        '{ op:assign_d_zero loc:<%intAdd(numStateVars,index)%> }'
    ;separator="\n")
    let assign_zero += (vars.algVars  |> var =>
        '{ op:assign_d_zero loc:<%intAdd(intMul(2,numStateVars),index)%> }'
    ;separator="\n")
    // states are independent variables
    let assign_ind = (vars.stateVars  |> var =>
        '{ op:assign_ind loc:<%index%> }'
    ;separator="\n")
    // derivates are dependent variables
    let assign_dep = (vars.derivativeVars  |> var =>
        '{ op:assign_dep loc:<%intAdd(numStateVars,index)%> }'
    ;separator="\n")
    
    let()= System.tmpTickResetIndex(0,28) /* reset ind index */
    let()= System.tmpTickIndexReserve(28, System.tmpTickMaximum(25))
    
    let equations = createEquation(odeEquations)
    
    let death_not = '{ op:death_not loc:0 loc:<%System.tmpTickMaximum(28)%> }'
    
    let text = assign_zero + assign_ind + equations + assign_dep + death_not;
    let()= textFile(text, '<%fileNamePrefix%>_adolcAsciiTrace.txt')
end generateAdolcAsciiTrace;

template createEquations(list<list<SimEqSystem>> eqs)
" create adolc equations operations"
::=
    let eqnops = (List.flatten(eqs) |> eq => 
        handle_equation(eq)
        ;separator="\n")
    << 
     <%eqnops%>
    >>
end createEquations;

template handle_equation(SimEqSystem eq)
"handle single equations for adolc trace"
::=
  match eq
  case SES_SIMPLE_ASSIGN(__) then
    let tmpOps = buffer ""
    let lhs = crefLoc(cref, tmpOps)
    let rhs = expLoc(exp, tmpOps)
end handle_equation;


template crefLoc(ComponentRef cref, Text &tmpOps)
"get cref location for adolc cref"
::=
  match cref
    case CREF_IDENT(ident = "time") then
        let timeLoc = System.tmpTickIndex(28)
        let tmpOps += '{ op:assign_p loc:0 loc:<%timeLoc%> }'
       'loc:<%timeLoc%>'
    else match cref2simvar(cr, getSimCode())
      case SIMVAR(varKind=PARAM()) then
        "(adouble)getparam(ad" + cref(cr) + ")"
      else
        "ad" + cref(cr)
end crefLoc;

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2