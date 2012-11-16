/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-CurrentYear, Linköping University,
 * Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 
 * AND THIS OSMC PUBLIC LICENSE (OSMC-PL). 
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES RECIPIENT'S  
 * ACCEPTANCE OF THE OSMC PUBLIC LICENSE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from Linköping University, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or  
 * http://www.openmodelica.org, and in the OpenModelica distribution. 
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS
 * OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

encapsulated package SCodeFlattenRedeclare
" file:        SCodeFlattenRedeclare.mo
  package:     SCodeFlattenRedeclare
  description: SCode flattening

  RCS: $Id$

  This module contains redeclare-specific functions used by SCodeFlatten to
  handle redeclares. There are three different types of redeclares that are
  handled: redeclare modifiers, element redeclares and class extends.

  REDECLARE MODIFIERS:
  Redeclare modifiers are redeclarations given as modifiers on an extends
  clause. When an extends clause is added to the environment with
  SCodeEnv.extendEnvWithExtends these modifiers are extracted with
  extractRedeclareFromModifier as a list of elements, and then stored in the
  SCodeEnv.EXTENDS representation. When SCodeLookup.lookupInBaseClasses is used
  to search for an identifier in a base class, these elements are replaced in
  the environment prior to searching in it by the replaceRedeclares function.

  ELEMENT REDECLARES:
  Element redeclares are similar to redeclare modifiers, but they are declared
  as standalone elements that redeclare an inherited element. When the
  environment is built they are initially added to a list of elements in the
  extends tables by SCodeEnv.addElementRedeclarationToEnvExtendsTable. When the
  environment is complete and SCodeEnv.updateExtendsInEnv is used to update the
  extends these redeclares are handled by addElementRedeclarationsToEnv, which
  looks up which base class each redeclare is redeclaring in. The element
  redeclares are then added to the list of redeclarations in the correct
  SCodeEnv.EXTENDS, and handled in the same way as redeclare modifiers.

  CLASS EXTENDS:
  Class extends are handled by adding them to the environment with
  extendEnvWithClassExtends. This function adds the given class as a normal
  class to the environment, and sets the class extends information field in
  the class's environment. This information is the base class and modifiers of
  the class extends. This information is later used when extends are updated
  with SCodeEnv.updateExtendsInEnv, and updateClassExtends is called.
  updateClassExtends looks up the full path to the base class of the class
  extends, and adds an extends clause to the class that extends from the base
  class. 
  
  However, since it's possible to redeclare the base class of a class
  extends it's possible that the base class is replaced with a class that
  extends from it. If the base class were to be replaced with this class it
  would mean that the class extends itself, causing a loop. To avoid this an
  alias for the base class is added instead, and the base class itself is added
  with the BASE_CLASS_SUFFIX defined in SCodeEnv. The alias can then be safely
  redeclared while preserving the base class for the class extends to extend
  from. It's somewhat difficult to only add aliases for classes that are used by
  class extends though, so an alias is added for all replaceable classes in
  SCodeEnv.extendEnvWithClassDef for simplicity's sake. The function
  SCodeLookup.resolveAlias is then used to resolve any alias items to the real
  items whenever an item is looked up in the environment.
  
  Class extends on the form 'redeclare class extends X' are thus
  translated to 'class X extends BaseClass.X$base', and then mostly handled like a
  normal class. Some care is needed in the dependency analysis to make sure
  that nothing important is removed, see comment in
  SCodeDependency.analyseClassExtends.  
"

public import Absyn;
public import SCode;
public import SCodeEnv;
public import InstTypes;
public import SCodeLookup;

public type Env = SCodeEnv.Env;
public type Item = SCodeEnv.Item;
public type Extends = SCodeEnv.Extends;
public type Prefix = InstTypes.Prefix;

public uniontype Replacement
  record REPLACED "an item got replaced"
    SCode.Ident name;
    Item old;
    Item new;
    Env env;
  end REPLACED;
  
  record PUSHED "the redeclares got pushed into the extends of the base classes"
    SCode.Ident name;
    Item redeclaredItem;
    list<Absyn.Path> baseClasses;
    SCodeEnv.ExtendsTable old;
    SCodeEnv.ExtendsTable new;
    Env env;
  end PUSHED;  
end Replacement;

public type Replacements = list<Replacement>;
public constant Replacements emptyReplacements = {};

protected import Debug;
protected import Error;
protected import Flags;
protected import List;
protected import SCodeCheck;
protected import Util;
protected import SCodeDump;

public function extendEnvWithClassExtends
  input SCode.Element inClassExtends;
  input Env inEnv;
  output Env outEnv;
algorithm
  outEnv := match(inClassExtends, inEnv)
    local
      SCode.Ident bc;
      list<SCode.Element> el;
      SCode.Partial pp;
      SCode.Encapsulated ep;
      SCode.Restriction res;
      SCode.Prefixes prefixes;
      Absyn.Info info;
      Env env, cls_env;
      SCode.Mod mods;
      Option<SCode.ExternalDecl> ext_decl;
      list<SCode.Annotation> annl;
      Option<SCode.Comment> cmt;
      list<SCode.Equation> nel, iel;
      list<SCode.AlgorithmSection> nal, ial;
      list<SCode.ConstraintSection> nco;
      list<Absyn.NamedArg> clats;  //class attributes
      SCode.ClassDef cdef;
      SCode.Element cls, ext;
      String el_str, env_str, err_msg;

    // When a 'redeclare class extends X' is encountered we insert a 'class X
    // extends BaseClass.X' into the environment, with the same elements as the
    // class extends clause. BaseClass is the class that class X is inherited
    // from. This allows us to look up elements in class extends, because
    // lookup can handle normal extends. This is the first phase where the
    // CLASS_EXTENDS is converted to a PARTS and added to the environment, and
    // the extends is added to the class environment's extends table. The
    // proper base class will be looked up in the second phase, in
    // updateClassExtends
    case (SCode.CLASS(
        prefixes = prefixes,
        encapsulatedPrefix = ep,
        partialPrefix = pp,
        restriction = res,
        classDef = SCode.CLASS_EXTENDS(
          baseClassName = bc, 
          modifications = mods,
          composition = SCode.PARTS(
            elementLst = el,
            normalEquationLst = nel,
            initialEquationLst = iel,
            normalAlgorithmLst = nal,
            initialAlgorithmLst = ial,
            constraintLst =  nco,
            clsattrs = clats,
            externalDecl = ext_decl,
            annotationLst = annl,
            comment = cmt)),
        info = info), _)
      equation
        // Construct a PARTS from the CLASS_EXTENDS.
        cdef = SCode.PARTS(el, nel, iel, nal, ial, nco, clats, ext_decl, annl, cmt);
        cls = SCode.CLASS(bc, prefixes, ep, pp, res, cdef, info);

        // Construct the class environment and add the new extends to it.
        cls_env = SCodeEnv.makeClassEnvironment(cls, false);
        ext = SCode.EXTENDS(Absyn.IDENT(bc), SCode.PUBLIC(), mods, NONE(), info);
        cls_env = addClassExtendsInfoToEnv(ext, cls_env);

        // Finally add the class to the environment.
        env = SCodeEnv.extendEnvWithItem(
          SCodeEnv.newClassItem(cls, cls_env, SCodeEnv.CLASS_EXTENDS()), inEnv, bc);
      then env;

    case (_, _)
      equation
        info = SCode.elementInfo(inClassExtends);
        el_str = SCodeDump.printElementStr(inClassExtends);
        env_str = SCodeEnv.getEnvName(inEnv);
        err_msg = "SCodeFlattenRedeclare.extendEnvWithClassExtends failed on unknown element " +& 
          el_str +& " in " +& env_str;
        Error.addSourceMessage(Error.INTERNAL_ERROR, {err_msg}, info);
      then
        fail();

  end match;
end extendEnvWithClassExtends;
  
protected function addClassExtendsInfoToEnv
  "Adds a class extends to the environment."
  input SCode.Element inClassExtends;
  input Env inEnv;
  output Env outEnv;
algorithm
  outEnv := matchcontinue(inClassExtends, inEnv)
    local
      list<Extends> bcl;
      list<SCode.Element> re;
      String estr;
      SCodeEnv.ExtendsTable ext;

    case (_, _)
      equation
        SCodeEnv.EXTENDS_TABLE(bcl, re, NONE()) = 
          SCodeEnv.getEnvExtendsTable(inEnv);
        ext = SCodeEnv.EXTENDS_TABLE(bcl, re, SOME(inClassExtends));
      then
        SCodeEnv.setEnvExtendsTable(ext, inEnv);

    else
      equation
        estr = "- SCodeFlattenRedeclare.addClassExtendsInfoToEnv: Trying to overwrite " +& 
               "existing class extends information, this should not happen!.";
        Error.addMessage(Error.INTERNAL_ERROR, {estr});
      then
        fail();

  end matchcontinue;
end addClassExtendsInfoToEnv;

public function updateClassExtends
  input SCode.Element inClass;
  input Env inEnv;
  input SCodeEnv.ClassType inClassType;
  output SCode.Element outClass;
  output Env outEnv;
algorithm
  (outClass, outEnv) := match(inClass, inEnv, inClassType)
    local
      String name;
      Env env;
      SCode.Mod mods;
      Absyn.Info info;
      SCode.Element cls, ext;

    case (_, SCodeEnv.FRAME(name = SOME(name), 
        extendsTable = SCodeEnv.EXTENDS_TABLE(classExtendsInfo = SOME(ext))) :: _,
        SCodeEnv.CLASS_EXTENDS())
      equation
        SCode.EXTENDS(modifications = mods, info = info) = ext;
        (cls, env) = updateClassExtends2(inClass, name, mods, info, inEnv);
      then
        (cls, env);

    else (inClass, inEnv);
  end match;
end updateClassExtends;

protected function updateClassExtends2
  input SCode.Element inClass;
  input String inName;
  input SCode.Mod inMods;
  input Absyn.Info inInfo;
  input Env inEnv;
  output SCode.Element outClass;
  output Env outEnv;
algorithm
  (outClass, outEnv) := match(inClass, inName, inMods, inInfo, inEnv)
    local
      Absyn.Path path;
      SCode.Element ext;
      SCodeEnv.Frame cls_frame;
      Env env;
      SCode.Element cls;
      Item item;

    case (_, _, _, _, cls_frame :: env)
      equation
        (path, item) = lookupClassExtendsBaseClass(inName, env, inInfo);
        SCodeCheck.checkClassExtendsReplaceability(item, Absyn.dummyInfo);
        ext = SCode.EXTENDS(path, SCode.PUBLIC(), inMods, NONE(), inInfo);
        {cls_frame} = SCodeEnv.extendEnvWithExtends(ext, {cls_frame});
        cls = SCode.addElementToClass(ext, inClass);
      then
        (cls, cls_frame :: env);

  end match;
end updateClassExtends2;

protected function lookupClassExtendsBaseClass
  "This function takes the name of a base class and looks up that name suffixed
   with the base class suffix defined in SCodeEnv. I.e. it looks up the real base
   class of a class extends, and not the alias introduced when adding replaceable
   classes to the environment in SCodeEnv.extendEnvWithClassDef. It returns the
   fully qualified path and the item for that base class."
  input String inName;
  input Env inEnv;
  input Absyn.Info inInfo;
  output Absyn.Path outPath;
  output Item outItem;
algorithm
  (outPath, outItem) := matchcontinue(inName, inEnv, inInfo)
    local
      Absyn.Path path;
      Item item;
      String basename;
      Env env;

    // Add the base class suffix to the name and try to look it up.
    case (_, _, _)
      equation
        basename = inName +& SCodeEnv.BASE_CLASS_SUFFIX;
        (item, _, env) = SCodeLookup.lookupInheritedName(basename, inEnv);
        // Use a special $ce qualified so that we can find the correct class
        // with SCodeLookup.lookupBaseClassName.
        path = Absyn.QUALIFIED("$ce", Absyn.IDENT(basename));
      then
        (path, item);

    // The previous case will fail if we try to class extend a
    // non-replaceable class, because they don't have aliases. To get the
    // correct error message later we look the class up via the non-alias name
    // instead and return that result if found.
    case (_, _, _)
      equation
        (item, path, env) = SCodeLookup.lookupInheritedName(inName, inEnv);
      then
        (path, item);
        
    else
      equation
        Error.addSourceMessage(Error.INVALID_REDECLARATION_OF_CLASS,
          {inName}, inInfo);
      then
        fail();

  end matchcontinue;
end lookupClassExtendsBaseClass;

public function addElementRedeclarationsToEnv
  input list<SCode.Element> inRedeclares;
  input Env inEnv;
  output Env outEnv;
algorithm
  outEnv := List.fold(inRedeclares, addElementRedeclarationsToEnv2, inEnv);
end addElementRedeclarationsToEnv;

protected function addElementRedeclarationsToEnv2
  input SCode.Element inRedeclare;
  input Env inEnv;
  output Env outEnv;
algorithm
  outEnv := matchcontinue(inRedeclare, inEnv)
    local
      SCode.Ident cls_name, name;
      Absyn.Info info;
      Absyn.Path ext_path, env_path;
      Env env;
      Item base_item, item;
      SCode.Element redecl;

    case (_, _)
      equation
        name = SCode.elementName(inRedeclare);
        info = SCode.elementInfo(inRedeclare);
        ext_path = lookupElementRedeclaration(name, inEnv, info);
        env_path = SCodeEnv.getEnvPath(inEnv);
        item = SCodeEnv.ALIAS(name, SOME(env_path), info);
        env = addRedeclareToEnvExtendsTable(item, ext_path, inEnv, info);
      then
        env;

    else
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.traceln("- SCodeFlattenRedeclare.addElementRedeclarationsToEnv failed for " +&
          SCode.elementName(inRedeclare) +& " in " +& 
          SCodeEnv.getEnvName(inEnv) +& "\n");
      then
        fail();
  end matchcontinue;
end addElementRedeclarationsToEnv2;

protected function lookupElementRedeclaration
  input SCode.Ident inName;
  input Env inEnv;
  input Absyn.Info inInfo;
  output Absyn.Path outPath;
algorithm
  outPath := matchcontinue(inName, inEnv, inInfo)
    local
      Absyn.Path path;

    case (_, _, _)
      equation
        path = SCodeLookup.lookupBaseClass(inName, inEnv);
      then
        path;

    else
      equation
        Error.addSourceMessage(Error.REDECLARE_NONEXISTING_ELEMENT,
          {inName}, inInfo);
      then
        fail();

  end matchcontinue;
end lookupElementRedeclaration;

protected function addRedeclareToEnvExtendsTable
  input Item inRedeclaredElement;
  input Absyn.Path inBaseClass;
  input Env inEnv;
  input Absyn.Info inInfo;
  output Env outEnv;
protected
  list<Extends> bcl;
  list<SCode.Element> re;
  Option<SCode.Element> cei;
algorithm
  SCodeEnv.EXTENDS_TABLE(bcl, re, cei) := SCodeEnv.getEnvExtendsTable(inEnv);
  bcl := addRedeclareToEnvExtendsTable2(inRedeclaredElement, inBaseClass, bcl);
  outEnv := SCodeEnv.setEnvExtendsTable(SCodeEnv.EXTENDS_TABLE(bcl, re, cei), inEnv);
end addRedeclareToEnvExtendsTable;

protected function addRedeclareToEnvExtendsTable2
  input Item inRedeclaredElement;
  input Absyn.Path inBaseClass;
  input list<Extends> inExtends;
  output list<Extends> outExtends;
algorithm
  outExtends := matchcontinue(inRedeclaredElement, inBaseClass, inExtends)
    local
      Extends ex;
      list<Extends> exl;
      Absyn.Path bc;
      list<SCodeEnv.Redeclaration> el;
      Absyn.Info info;
      SCodeEnv.Redeclaration redecl;

    case (_, _, SCodeEnv.EXTENDS(bc, el, info) :: exl)
      equation
        true = Absyn.pathEqual(inBaseClass, bc);
        redecl = SCodeEnv.PROCESSED_MODIFIER(inRedeclaredElement);
        SCodeCheck.checkDuplicateRedeclarations(redecl, el);
        ex = SCodeEnv.EXTENDS(bc, redecl :: el, info);
      then
        ex :: exl;

    case (_, _, ex :: exl)
      equation
        exl = addRedeclareToEnvExtendsTable2(inRedeclaredElement, inBaseClass, exl);
      then
        ex :: exl;
    
  end matchcontinue;
end addRedeclareToEnvExtendsTable2;

public function processRedeclare
  "Processes a raw redeclare modifier into a processed form."
  input SCodeEnv.Redeclaration inRedeclare;
  input Env inEnv;
  input InstTypes.Prefix inPrefix;
  output SCodeEnv.Redeclaration outRedeclare;
algorithm
  outRedeclare := matchcontinue(inRedeclare, inEnv, inPrefix)
    local
      SCode.Ident name;
      SCode.Partial pp;
      SCode.Encapsulated ep;
      SCode.Prefixes prefixes;
      Absyn.Path path;
      SCode.Mod mod;
      Option<SCode.Comment> cmt;
      SCode.Restriction res;
      Absyn.Info info;
      SCode.Attributes attr;
      Option<Absyn.Exp> cond;
      Option<Absyn.ArrayDim> ad;

      Item el_item, redecl_item;
      SCode.Element el;
      Env cls_env, env;
   
   case (SCodeEnv.RAW_MODIFIER(modifier = el as SCode.CLASS(name = _)), _, _)
      equation
        cls_env = SCodeEnv.makeClassEnvironment(el, true);
        el_item = SCodeEnv.newClassItem(el, cls_env, SCodeEnv.USERDEFINED());
        redecl_item = SCodeEnv.REDECLARED_ITEM(el_item, inEnv);
      then
        SCodeEnv.PROCESSED_MODIFIER(redecl_item);

    case (SCodeEnv.RAW_MODIFIER(modifier = el as SCode.COMPONENT(name = _)), _, _)
      equation
        el_item = SCodeEnv.newVarItem(el, true);
        redecl_item = SCodeEnv.REDECLARED_ITEM(el_item, inEnv);
      then
        SCodeEnv.PROCESSED_MODIFIER(redecl_item);

    case (SCodeEnv.PROCESSED_MODIFIER(modifier = _), _, _) then inRedeclare;

    else
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.traceln("- SCodeFlattenRedeclare.processRedeclare failed on " +&
          SCodeDump.printElementStr(SCodeEnv.getRedeclarationElement(inRedeclare)) +& 
          " in " +& Absyn.pathString(SCodeEnv.getEnvPath(inEnv)));
      then
        fail();
  end matchcontinue;
end processRedeclare;

public function replaceRedeclares
  "Replaces redeclares in the environment. This function takes a list of
   redeclares, the item and environment of the class in which they should be
   redeclared, and the environment in which the modified element was declared
   (used to qualify the redeclares). The redeclares are then either replaced if
   they can be found in the immediate local environment of the class, or pushed
   into the correct extends clauses if they are inherited." 
  input list<SCodeEnv.Redeclaration> inRedeclares;
  input Item inClassItem "The item of the class to be modified.";
  input Env inClassEnv "The environment of the class to be modified.";
  input Env inElementEnv "The environment in which the modified element was declared.";
  input SCodeLookup.RedeclareReplaceStrategy inReplaceRedeclares;
  output Option<Item> outItem;
  output Option<Env> outEnv;
algorithm
  (outItem, outEnv) := matchcontinue(inRedeclares, inClassItem, inClassEnv,
      inElementEnv, inReplaceRedeclares)
    local
      Item item;
      Env env;

    case (_, _, _, _, SCodeLookup.IGNORE_REDECLARES()) 
      then (SOME(inClassItem), SOME(inClassEnv));

    case (_, _, _, _, SCodeLookup.INSERT_REDECLARES())
      equation
        (item, env, _) = replaceRedeclaredElementsInEnv(inRedeclares,
          inClassItem, inClassEnv, inElementEnv, InstTypes.emptyPrefix);
      then
        (SOME(item), SOME(env));

    else (NONE(), NONE());
  end matchcontinue;
end replaceRedeclares;

public function replaceRedeclaredElementsInEnv
  "If a variable or extends clause has modifications that redeclare classes in
   it's instance we need to replace those classes in the environment so that the
   lookup finds the right classes. This function takes a list of redeclares from
   an elements' modifications and applies them to the environment of the
   elements type."
  input list<SCodeEnv.Redeclaration> inRedeclares "The redeclares from the modifications.";
  input Item inItem "The type of the element.";
  input Env inTypeEnv "The enclosing scopes of the type.";
  input Env inElementEnv "The environment in which the element was declared.";
  input InstTypes.Prefix inPrefix;
  output Item outItem;
  output Env outEnv;
  output Replacements outReplacements "what replacements where performed if any";
algorithm
  (outItem, outEnv, outReplacements) :=
  matchcontinue(inRedeclares, inItem, inTypeEnv, inElementEnv, inPrefix)
    local
      SCode.Element cls;
      Env env;
      SCodeEnv.Frame item_env;
      SCodeEnv.ClassType cls_ty;
      list<SCodeEnv.Redeclaration> redecls;
      Replacements repl;

    // no redeclares!
    case ({}, _, _, _, _) then (inItem, inTypeEnv, {});

    case (_, SCodeEnv.CLASS(cls = cls, env = {item_env}, classType = cls_ty), _, _, _)
      equation
        // Merge the types environment with it's enclosing scopes to get the
        // enclosing scopes of the classes we need to replace.
        env = SCodeEnv.enterFrame(item_env, inTypeEnv);
        redecls = List.map2(inRedeclares, processRedeclare, inElementEnv, inPrefix);
        ((env, repl)) = List.fold(redecls, replaceRedeclaredElementInEnv, ((env, emptyReplacements)));
        item_env :: env = env;
      then
        (SCodeEnv.CLASS(cls, {item_env}, cls_ty), env, repl);

    else
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.trace("- SCodeFlattenRedeclare.replaceRedeclaredElementsInEnv failed for: ");
        Debug.traceln("redeclares: " +& 
          stringDelimitList(List.map(inRedeclares, SCodeEnv.printRedeclarationStr), "\n---------\n") +&  
          " item: " +& SCodeEnv.getItemName(inItem) +& " in scope:" +& SCodeEnv.getEnvName(inElementEnv));
      then
        fail();
  end matchcontinue;
end replaceRedeclaredElementsInEnv;

public function extractRedeclaresFromModifier
  "Returns a list of redeclare elements given a redeclaration modifier."
  input SCode.Mod inMod;
  output list<SCodeEnv.Redeclaration> outRedeclares;
algorithm
  outRedeclares := match(inMod)
    local
      list<SCode.SubMod> sub_mods;
      list<SCodeEnv.Redeclaration> redeclares;
    
    case SCode.MOD(subModLst = sub_mods)
      equation
        redeclares = List.fold(sub_mods, extractRedeclareFromSubMod, {});
      then
        redeclares;

    else then {};
  end match;
end extractRedeclaresFromModifier;

protected function extractRedeclareFromSubMod
  "Checks a submodifier and adds the redeclare element to the list of redeclares
  if the modifier is a redeclaration modifier."
  input SCode.SubMod inMod;
  input list<SCodeEnv.Redeclaration> inRedeclares;
  output list<SCodeEnv.Redeclaration> outRedeclares;
algorithm
  outRedeclares := match(inMod, inRedeclares)
    local
      SCode.Element el;
      SCodeEnv.Redeclaration redecl; 

    case (SCode.NAMEMOD(A = SCode.REDECL(element = el)), _)
      equation
        redecl = SCodeEnv.RAW_MODIFIER(el);
        SCodeCheck.checkDuplicateRedeclarations(redecl, inRedeclares);
      then
        redecl :: inRedeclares;

    // Skip modifiers that are not redeclarations.
    else inRedeclares;
  end match;
end extractRedeclareFromSubMod;

protected function replaceRedeclaredElementInEnv
  "Replaces a redeclaration in the environment."
  input SCodeEnv.Redeclaration inRedeclare;
  input tuple<Env, Replacements> inEnv;
  output tuple<Env, Replacements> outEnv;
algorithm
  outEnv := matchcontinue(inRedeclare, inEnv)
    local
      SCode.Ident name, scope_name;
      Item item;
      Absyn.Info info;
      list<Absyn.Path> bcl;
      list<String> bcl_str;
      Env env;
      tuple<Env, Replacements> envRpl;

    // Try to redeclare this element in the current scope.
    case (SCodeEnv.PROCESSED_MODIFIER(modifier = item), _)
      equation
        name = SCodeEnv.getItemName(item);
        // do not asume the story ends here
        // you have to push into extends again 
        // even if you find it in the local scope!
        envRpl = pushRedeclareIntoExtendsNoFail(name, item, inEnv);
      then  
        replaceElementInScope(name, item, envRpl);
        
    // If the previous case failed, see if we can find the redeclared element in
    // any of the base classes. If so, push the redeclare into those base
    // classes instead, i.e. add them to the list of redeclares in the
    // appropriate extends in the extends table.
    case (SCodeEnv.PROCESSED_MODIFIER(modifier = item), _)
      equation
        name = SCodeEnv.getItemName(item);
        bcl = SCodeLookup.lookupBaseClasses(name, Util.tuple21(inEnv));
      then
        pushRedeclareIntoExtends(name, item, bcl, inEnv);
        
    // The redeclared element could not be found, show an error.
    case (SCodeEnv.PROCESSED_MODIFIER(modifier = item), _)
      equation
        scope_name = SCodeEnv.getScopeName(Util.tuple21(inEnv));
        name = SCodeEnv.getItemName(item);
        info = SCodeEnv.getItemInfo(item);
        Error.addSourceMessage(Error.MISSING_MODIFIED_ELEMENT,
          {name, scope_name}, info);
      then
        fail(); 

  end matchcontinue;
end replaceRedeclaredElementInEnv;

protected function pushRedeclareIntoExtendsNoFail
"Pushes a redeclare into the given extends in the environment if it can.
 if not just returns the same tuple<env, repl>"
  input SCode.Ident inName;
  input Item inRedeclare;
  input tuple<Env, Replacements> inEnv;
  output tuple<Env, Replacements> outEnv;
algorithm
  outEnv := matchcontinue(inName, inRedeclare, inEnv)
    local
      SCode.Ident name, scope_name;
      Item item;
      Absyn.Info info;
      list<Absyn.Path> bcl;
      list<String> bcl_str;
      Env env;
      tuple<Env, Replacements> envRpl;
    
    case (_, _, _)
      equation
        bcl = SCodeLookup.lookupBaseClasses(inName, Util.tuple21(inEnv));
        (envRpl) = pushRedeclareIntoExtends(inName, inRedeclare, bcl, inEnv);
      then
        envRpl;
    
    else inEnv;
  end matchcontinue;
end pushRedeclareIntoExtendsNoFail;

protected function pushRedeclareIntoExtends
  "Pushes a redeclare into the given extends in the environment."
  input SCode.Ident inName;
  input Item inRedeclare;
  input list<Absyn.Path> inBaseClasses;
  input tuple<Env, Replacements> inEnv;
  output tuple<Env, Replacements> outEnv;
protected
  list<SCodeEnv.Extends> exts;
  list<SCode.Element> re;
  Option<SCode.Element> cei;
  SCodeEnv.ExtendsTable etNew, etOld;
  String name;
  Env env;
  Replacements repl;
algorithm
  (env, repl) := inEnv; 
  
  SCodeEnv.FRAME(extendsTable = etOld as SCodeEnv.EXTENDS_TABLE(exts, re, cei)) :: _ := env;
  exts := pushRedeclareIntoExtends2(inName, inRedeclare, inBaseClasses, exts);
  etNew := SCodeEnv.EXTENDS_TABLE(exts, re, cei);
  
  env := SCodeEnv.setEnvExtendsTable(etNew, env);
  repl := PUSHED(inName, inRedeclare, inBaseClasses, etOld, etNew, env)::repl;
  
  outEnv := (env, repl);
  // tracePushRedeclareIntoExtends(inName, inRedeclare, inBaseClasses, env, etOld, etNew);
end pushRedeclareIntoExtends;

protected function pushRedeclareIntoExtends2
  "This function takes a redeclare item and a list of base class paths that the
   redeclare item should be added to. It goes through the given list of
   extends and pushes the redeclare into each one that's in the list of the
   base class paths. It assumes that the list of base class paths and extends
   are sorted in the same order."
  input String inName;
  input Item inRedeclare;
  input list<Absyn.Path> inBaseClasses;
  input list<SCodeEnv.Extends> inExtends;
  output list<SCodeEnv.Extends> outExtends;
algorithm
  outExtends := matchcontinue(inName, inRedeclare, inBaseClasses, inExtends)
    local
      Absyn.Path bc1, bc2;
      list<Absyn.Path> rest_bc;
      SCodeEnv.Extends ext;
      list<SCodeEnv.Extends> rest_exts;
      list<SCodeEnv.Redeclaration> redecls;
      Absyn.Info info;
      list<String> bc_strl;
      String bcl_str, err_msg;

    // See if the first base class path matches the first extends. Push the
    // redeclare into that extends if so.
    case (_, _, bc1 :: rest_bc, SCodeEnv.EXTENDS(bc2, redecls, info) :: rest_exts)
      equation
        true = Absyn.pathEqual(bc1, bc2);
        redecls = pushRedeclareIntoExtends3(inRedeclare, inName, redecls);
        rest_exts = pushRedeclareIntoExtends2(inName, inRedeclare, rest_bc, rest_exts);
      then
        SCodeEnv.EXTENDS(bc2, redecls, info) :: rest_exts;

    // The extends didn't match, continue with the rest of them.
    case (_, _, rest_bc, ext :: rest_exts)
      equation
        rest_exts = pushRedeclareIntoExtends2(inName, inRedeclare, rest_bc, rest_exts);
      then
        ext :: rest_exts;

    // No more base class paths to match means we're done.
    case (_, _, {}, _) then inExtends;

    // No more extends means that we couldn't find all the base classes. This
    // shouldn't happen.
    case (_, _, _, {})
      equation
        bc_strl = List.map(inBaseClasses, Absyn.pathString);
        bcl_str = stringDelimitList(bc_strl, ", ");
        err_msg = "SCodeFlattenRedeclare.pushRedeclareIntoExtends2 couldn't find the base classes {"
          +& bcl_str +& "} for " +& inName;
        Error.addMessage(Error.INTERNAL_ERROR, {err_msg});
      then
        fail();

  end matchcontinue;
end pushRedeclareIntoExtends2;

protected function pushRedeclareIntoExtends3
  "Given the item and name of a redeclare, try to find the redeclare in the
   given list of redeclares. If found, replace the redeclare in the list.
   Otherwise, add a new redeclare to the list."
  input Item inRedeclare;
  input String inName;
  input list<SCodeEnv.Redeclaration> inRedeclares;
  output list<SCodeEnv.Redeclaration> outRedeclares;
algorithm
  outRedeclares := matchcontinue(inRedeclare, inName, inRedeclares)
    local
      Item item;
      SCodeEnv.Redeclaration redecl;
      list<SCodeEnv.Redeclaration> rest_redecls;
      String name;

    case (_, _, SCodeEnv.PROCESSED_MODIFIER(modifier = item) :: rest_redecls)
      equation
        name = SCodeEnv.getItemName(item);
        true = stringEqual(name, inName);
      then
        SCodeEnv.PROCESSED_MODIFIER(inRedeclare) :: rest_redecls;

    case (_, _, redecl :: rest_redecls)
      equation
        rest_redecls = pushRedeclareIntoExtends3(inRedeclare, inName, rest_redecls);
      then
        redecl :: rest_redecls;

    case (_, _, {}) then {SCodeEnv.PROCESSED_MODIFIER(inRedeclare)};

  end matchcontinue;
end pushRedeclareIntoExtends3;
        
public function replaceElementInScope
  "Replaces an element in the current scope."
  input SCode.Ident inElementName;
  input Item inElement;
  input tuple<Env, Replacements> inEnv;
  output tuple<Env, Replacements> outEnv;
algorithm
  outEnv := match(inElementName, inElement, inEnv)
    local
      SCodeEnv.AvlTree tree;
      Item old_item, new_item;
      Env env;
      Replacements repl;

    case (_, _, (env as SCodeEnv.FRAME(clsAndVars = tree) :: _, repl))
      equation
        old_item = SCodeEnv.avlTreeGet(tree, inElementName);
        /*********************************************************************/
        // TODO: Check if this is actually needed
        /*********************************************************************/
        new_item = propagateItemPrefixes(old_item, inElement);
        new_item = SCodeEnv.linkItemUsage(old_item, new_item);
        tree = SCodeEnv.avlTreeReplace(tree, inElementName, new_item);
        env = SCodeEnv.setEnvClsAndVars(tree, env);
        repl = REPLACED(inElementName, old_item, new_item, env)::repl;
        // traceReplaceElementInScope(inElementName, old_item, new_item, env);
      then
        ((env, repl));

  end match;
end replaceElementInScope;

protected function propagateItemPrefixes
  input Item inOriginalItem;
  input Item inNewItem;
  output Item outNewItem;
algorithm
  outNewItem := match(inOriginalItem, inNewItem)
    local
      SCode.Element el1, el2;
      Option<Util.StatefulBoolean> iu1, iu2;
      Env env1, env2;
      SCodeEnv.ClassType ty1, ty2;
      Item item;

    case (SCodeEnv.VAR(var = el1, isUsed = iu1), 
          SCodeEnv.VAR(var = el2, isUsed = iu2))
      equation
        el2 = propagateAttributesVar(el1, el2);
      then
        SCodeEnv.VAR(el2, iu2);

    case (SCodeEnv.CLASS(cls = el1, env = env1, classType = ty1),
          SCodeEnv.CLASS(cls = el2, env = env2, classType = ty2))
      equation
        el2 = propagateAttributesClass(el1, el2);
      then
        SCodeEnv.CLASS(el2, env2, ty2);

    /*************************************************************************/
    // TODO: Attributes should probably be propagated for alias items too. If
    // the original is an alias, look up the referenced item and use those
    // attributes. If the new item is an alias, look up the referenced item and
    // apply the attributes to it.
    /*************************************************************************/
    case (SCodeEnv.ALIAS(path = _), _) then inNewItem;
    case (_, SCodeEnv.ALIAS(path = _)) then inNewItem;

    case (SCodeEnv.REDECLARED_ITEM(item = item), _)
      then propagateItemPrefixes(item, inNewItem);

    case (_, SCodeEnv.REDECLARED_ITEM(item = item, declaredEnv = env1))
      equation
        item = propagateItemPrefixes(inOriginalItem, item);
      then
      SCodeEnv.REDECLARED_ITEM(item, env1);

    else
      equation
        Error.addMessage(Error.INTERNAL_ERROR,
          {"SCodeFlattenRedeclare.propagateAttributes failed on unknown item."});
      then
        fail();
  end match;
end propagateItemPrefixes;

protected function propagateAttributesVar
  input SCode.Element inOriginalVar;
  input SCode.Element inNewVar;
  output SCode.Element outNewVar;
protected
  SCode.Ident name;
  SCode.Prefixes pref1, pref2;
  SCode.Attributes attr1, attr2;
  Absyn.TypeSpec ty;
  SCode.Mod mod;
  Option<SCode.Comment> cmt;
  Option<Absyn.Exp> cond;
  Absyn.Info info;
algorithm
  SCode.COMPONENT(prefixes = pref1, attributes = attr1) := inOriginalVar;
  SCode.COMPONENT(name, pref2, attr2, ty, mod, cmt, cond, info) := inNewVar;
  pref2 := propagatePrefixes(pref1, pref2);
  attr2 := propagateAttributes(attr1, attr2);
  outNewVar := SCode.COMPONENT(name, pref2, attr2, ty, mod, cmt, cond, info);
end propagateAttributesVar;

public function propagateAttributesClass
  input SCode.Element inOriginalClass;
  input SCode.Element inNewClass;
  output SCode.Element outNewClass;
protected
  SCode.Ident name;
  SCode.Prefixes pref1, pref2;
  SCode.Encapsulated ep;
  SCode.Partial pp;
  SCode.Restriction res;
  SCode.ClassDef cdef;
  Absyn.Info info;
algorithm
  SCode.CLASS(prefixes = pref1) := inOriginalClass;
  SCode.CLASS(name, pref2, ep, pp, res, cdef, info) := inNewClass;
  pref2 := propagatePrefixes(pref1, pref2);
  outNewClass := SCode.CLASS(name, pref2, ep, pp, res, cdef, info);
end propagateAttributesClass;
    
protected function propagatePrefixes
  input SCode.Prefixes inOriginalPrefixes;
  input SCode.Prefixes inNewPrefixes;
  output SCode.Prefixes outNewPrefixes;
protected
  SCode.Visibility vis1, vis2;
  Absyn.InnerOuter io1, io2;
  SCode.Redeclare rdp;
  SCode.Final fp;
  SCode.Replaceable rpp;
algorithm
  SCode.PREFIXES(visibility = vis1, innerOuter = io1) := inOriginalPrefixes;
  SCode.PREFIXES(vis2, rdp, fp, io2, rpp) := inNewPrefixes;
  io2 := propagatePrefixInnerOuter(io1, io2);
  outNewPrefixes := SCode.PREFIXES(vis2, rdp, fp, io2, rpp);
end propagatePrefixes;

protected function propagatePrefixInnerOuter
  input Absyn.InnerOuter inOriginalIO;
  input Absyn.InnerOuter inIO;
  output Absyn.InnerOuter outIO;
algorithm
  outIO := match(inOriginalIO, inIO)
    case (_, Absyn.NOT_INNER_OUTER()) then inOriginalIO;
    else inIO;
  end match;
end propagatePrefixInnerOuter;

protected function propagateAttributes
  input SCode.Attributes inOriginalAttributes;
  input SCode.Attributes inNewAttributes;
  output SCode.Attributes outNewAttributes;
protected
  Absyn.ArrayDim dims1, dims2;
  SCode.ConnectorType ct1, ct2;
  SCode.Parallelism prl1,prl2;
  SCode.Variability var1, var2;
  Absyn.Direction dir1, dir2;
algorithm
  SCode.ATTR(dims1, ct1, prl1, var1, dir1) := inOriginalAttributes;
  SCode.ATTR(dims2, ct2, prl2, var2, dir2) := inNewAttributes;
  dims2 := propagateArrayDimensions(dims1, dims2);
  ct2 := propagateConnectorType(ct1, ct2);
  prl2 := propagateParallelism(prl1,prl2);
  var2 := propagateVariability(var1, var2);
  dir2 := propagateDirection(dir1, dir2);
  outNewAttributes := SCode.ATTR(dims2, ct2, prl2, var2, dir2);
end propagateAttributes;

protected function propagateArrayDimensions
  input Absyn.ArrayDim inOriginalDims;
  input Absyn.ArrayDim inNewDims;
  output Absyn.ArrayDim outNewDims;
algorithm
  outNewDims := match(inOriginalDims, inNewDims)
    case (_, {}) then inOriginalDims;
    else inNewDims;
  end match;
end propagateArrayDimensions;

protected function propagateConnectorType
  input SCode.ConnectorType inOriginalConnectorType;
  input SCode.ConnectorType inNewConnectorType;
  output SCode.ConnectorType outNewConnectorType;
algorithm
  outNewConnectorType := match(inOriginalConnectorType, inNewConnectorType)
    case (_, SCode.POTENTIAL()) then inOriginalConnectorType;
    else inNewConnectorType;
  end match;
end propagateConnectorType;

protected function propagateParallelism
  input SCode.Parallelism inOriginalParallelism;
  input SCode.Parallelism inNewParallelism;
  output SCode.Parallelism outNewParallelism;
algorithm
  outNewParallelism := match(inOriginalParallelism, inNewParallelism)
    case (_, SCode.NON_PARALLEL()) then inOriginalParallelism;
    else inNewParallelism;
  end match;
end propagateParallelism;

protected function propagateVariability
  input SCode.Variability inOriginalVariability;
  input SCode.Variability inNewVariability;
  output SCode.Variability outNewVariability;
algorithm
  outNewVariability := match(inOriginalVariability, inNewVariability)
    case (_, SCode.VAR()) then inOriginalVariability;
    else inNewVariability;
  end match;
end propagateVariability;

protected function propagateDirection
  input Absyn.Direction inOriginalDirection;
  input Absyn.Direction inNewDirection;
  output Absyn.Direction outNewDirection;
algorithm
  outNewDirection := match(inOriginalDirection, inNewDirection)
    case (_, Absyn.BIDIR()) then inOriginalDirection;
    else inNewDirection;
  end match;
end propagateDirection;

protected function traceReplaceElementInScope
"@author: adrpo
 good for debugging redeclares.
 uncomment it in replaceElementInScope to activate it"
  input SCode.Ident inElementName;
  input Item inOldItem;
  input Item inNewItem;
  input Env inEnv;
algorithm
  _ := matchcontinue(inElementName, inOldItem, inNewItem, inEnv)
    case (_, _, _, _)
      equation
        print("replacing element: " +& inElementName +& " env: " +& SCodeEnv.getEnvName(inEnv) +& "\n\t");
        print("Old Element:" +& SCodeEnv.itemStr(inOldItem) +& 
              " env: " +& SCodeEnv.getEnvName(SCodeEnv.getItemEnvNoFail(inOldItem)) +& "\n\t");
        print("New Element:" +& SCodeEnv.itemStr(inNewItem) +& 
              " env: " +& SCodeEnv.getEnvName(SCodeEnv.getItemEnvNoFail(inNewItem)) +& 
              "\n===============\n");
      then ();
    
    else
      equation
        print("traceReplaceElementInScope failed on element: " +& inElementName +& "\n");
      then ();
  end matchcontinue;
end traceReplaceElementInScope;

protected function tracePushRedeclareIntoExtends
"@author: adrpo
 good for debugging redeclares.
 uncomment it in pushRedeclareIntoExtends to activate it"
  input SCode.Ident inName;
  input SCodeEnv.Item inRedeclare;
  input list<Absyn.Path> inBaseClasses;
  input Env inEnv;
  input SCodeEnv.ExtendsTable inEtNew;
  input SCodeEnv.ExtendsTable inEtOld;
algorithm
  _ := matchcontinue(inName, inRedeclare, inBaseClasses, inEnv, inEtNew, inEtOld)
    case (_, _, _, _, _, _)
      equation
        print("pushing: " +& inName +& " redeclare: " +& SCodeEnv.itemStr(inRedeclare) +& "\n\t"); 
        print("into baseclases: " +& stringDelimitList(List.map(inBaseClasses, Absyn.pathString), ", ") +& "\n\t");
        print("called from env: " +& SCodeEnv.getEnvName(inEnv) +& "\n");
        print("-----------------\n");
      then ();

    else
      equation
        print("tracePushRedeclareIntoExtends failed on element: " +& inName +& "\n");
      then ();

  end matchcontinue;
end tracePushRedeclareIntoExtends;

end SCodeFlattenRedeclare;
