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
  case simCode as SIMCODE(modelInfo=MODELINFO(vars=vars as SIMVARS(__),
                                              varInfo=varInfo as VARINFO(__)),
                          odeEquations=odeEquations) then
    let()= System.tmpTickResetIndex(0,25) /* reset tmp index */
    // states are independent variables
    let assign_zero = ""
    let &assign_zero += (vars.stateVars |> var as  SIMVAR(__) =>
            '{ op:assign_d_zero loc:<%index%> }'
    ;separator="\n")

    let &assign_zero += (vars.derivativeVars  |> var as SIMVAR(__) =>
        '{ op:assign_d_zero loc:<%intAdd(varInfo.numStateVars,index)%> }'
    ;separator="\n")
    let &assign_zero += (vars.algVars  |> var as SIMVAR(__) =>
        '{ op:assign_d_zero loc:<%intAdd(intMul(2,varInfo.numStateVars),index)%> }'
    ;separator="\n")
    // states are independent variables
    let assign_ind = ""
    let &assign_ind += (vars.stateVars  |> var as SIMVAR(__) =>
        '{ op:assign_ind loc:<%index%> }'
    ;separator="\n")
    // derivates are dependent variables
    let assign_dep = ""
    let &assign_dep += (vars.derivativeVars  |> var as SIMVAR(__) =>
        '{ op:assign_dep loc:<%intAdd(varInfo.numStateVars,index)%> }'
    ;separator="\n")
    
    let()= System.tmpTickResetIndex(0,28) /* reset ind index */
    let tickMax25 = System.tmpTickIndexReserve(28, System.tmpTickMaximum(25))
    
    let equations = createEquations(odeEquations)
    
    let death_not = '{ op:death_not loc:0 loc:<%System.tmpTickMaximum(28)%> }'
    
    let text = assign_zero + assign_ind + equations + assign_dep + death_not
    let()=textFile(text, '<%fileNamePrefix%>_adolcAsciiTrace.txt')
    <<
    // ADOLC trace file in ascii format
    // states are independent variables
    <%assign_zero%>
    // states are independent variables
    <%assign_ind%>
    // derivates are dependent variables
    <%assign_dep%>
    >>
  end match
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
    let tmpOps = ""
    let lhs = crefLoc(cref, tmpOps)
    let rhs = expLoc(exp, tmpOps)
  <<
  <%tmpOps%>
  <%lhs%> = <%rhs%>
  >>
  end match
end handle_equation;


template crefLoc(ComponentRef cref, Text &tmpOps)
"get cref location for adolc cref"
::=
  match cref
    case CREF_IDENT(ident = "time") then
      let timeLoc = System.tmpTickIndex(28)
      let &tmpOps += '{ op:assign_p loc:0 loc:<%timeLoc%> }'
      '<loc:<%timeLoc%>>'
    else match cref2simvar(cref, getSimCode())
      case SIMVAR(varKind=PARAM()) then
        "(adouble)getparam(ad" + crefStr(cref) + ")"
      else
        "ad" + crefStr(cref)
      end match
  end match
end crefLoc;

template expLoc(DAE.Exp exp, Text &tmpOps)
"get exp location for adolc cref"
::=
  ""
end expLoc;

annotation(__OpenModelica_Interface="backend");
end CodegenADOLC;

// vim: filetype=susan sw=2 sts=2