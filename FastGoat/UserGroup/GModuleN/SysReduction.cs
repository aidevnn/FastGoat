using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public record SysReduction(Polynomial<ZnInt, Xi> eq, 
    Polynomial<ZnInt, Xi> expr,
    int order,
    Xi xi, 
    Polynomial<ZnInt, Xi> substitutionExpr, 
    int mod);