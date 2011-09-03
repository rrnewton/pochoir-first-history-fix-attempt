{-
 ----------------------------------------------------------------------------------
 -  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 -  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 - 		                     Charles E. Leiserson <cel@mit.edu>
 - 	 
 -   This program is free software: you can redistribute it and/or modify
 -   it under the terms of the GNU General Public License as published by
 -   the Free Software Foundation, either version 3 of the License, or
 -   (at your option) any later version.
 -
 -   This program is distributed in the hope that it will be useful,
 -   but WITHOUT ANY WARRANTY; without even the implied warranty of
 -   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -   GNU General Public License for more details.
 -
 -   You should have received a copy of the GNU General Public License
 -   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 -
 -   Suggestsions:                  yuantang@csail.mit.edu
 -   Bugs:                          yuantang@csail.mit.edu
 -
 --------------------------------------------------------------------------------
 -}

-- The Basic Parser for a second pass --
module PBasicParser2 where

import Text.ParserCombinators.Parsec
import qualified Text.ParserCombinators.Parsec.Token as Token
import Text.ParserCombinators.Parsec.Expr
import Text.ParserCombinators.Parsec.Language
import Text.Read (read)

import Data.Char
import Data.List
import qualified Data.Map as Map

import PBasicParser
import PShow
import PUtils
import PData

ppStencil1 :: String -> ParserState -> GenParser Char ParserState String
ppStencil1 l_id l_state = 
    -- convert "Register_Kernel(g, k, ... ks) " to "Register_Obase_Kernel(g, k, bk)"
        do try $ pMember "Register_Stagger_Kernels"
           (l_guard, l_kernels) <- parens pStencilRegisterStaggerKernelParams
           semi
           case Map.lookup l_id $ pStencil l_state of
               Nothing -> return (l_id ++ ".Register_Stagger_Kernels(" ++ 
                                  (gName l_guard) ++ (gComment l_guard) ++ ", " ++ 
                                  (intercalate ", " $ zipWith (++) 
                                       (map kName l_kernels) 
                                       (map kComment l_kernels)) ++ 
                                  "); /* UNKNOWN Stencil " ++ l_id ++ " */" ++ 
                                  breakline)
               Just l_stencil ->
                   do let l_regKernels = pShowRegStaggerKernel (pMode l_state) l_stencil (l_guard, l_kernels)
                      return (l_regKernels)
    -- convert "Register_Kernel(g, k, ... ks) " to "Register_Obase_Kernel(g, k, bk)"
    <|> do try $ pMember "Register_Tile_Kernels"
           (l_guard, l_tile) <- parens pStencilRegisterTileKernelParams
           semi
           case Map.lookup l_id $ pStencil l_state of
               Nothing -> return (l_id ++ ".Register_Tile_Kernels(" ++ 
                                  (gName l_guard) ++ (gComment l_guard) ++ ", " ++
                                  (tName l_tile) ++ (tComment l_tile) ++
                                  "); /* UNKNOWN Stencil " ++ l_id ++ " */" ++ 
                                  breakline)
               Just l_stencil ->
               -- convert "Register_Kernel(g, k, ... ks) " to "Register_Obase_Kernel(g, k, bk)"
                   do let l_regKernels = pShowRegTileKernel (pMode l_state) l_stencil (l_guard, l_tile)
                      return (l_regKernels)
    -- Ad hoc implementation of Run_Unroll
    <|> do try $ pMember "Run"
           l_tstep <- parens exprStmtDim
           semi
           let l_mode = pMode l_state
           case Map.lookup l_id $ pStencil l_state of
               Nothing -> return (l_id ++ ".Run(" ++ show l_tstep ++ ");")
               Just l_stencil -> if l_mode == PMUnroll 
                                    || l_mode == PAllCondTileMacro 
                                    || l_mode == PAllCondTileCPointer 
                                    || l_mode == PAllCondTilePointer
                                    || l_mode == PAllCondTileOptPointer
                                    then return (breakline ++ l_id ++ 
                                         ".Run_Obase_Merge(" ++ show l_tstep ++ 
                                         "); /* Run with Stencil " ++ l_id ++ " */" ++ 
                                         breakline)
                                    else return (breakline ++ l_id ++ ".Run_Obase(" ++
                                         show l_tstep ++ "); /* Run with Stencil " ++
                                         l_id ++ " */" ++ breakline)
    <|> do return (l_id)

