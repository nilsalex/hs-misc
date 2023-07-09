module Main (main) where

import Control.Monad.Random.Lazy (MonadRandom (..), MonadTrans (lift), RandT, RandomGen, evalRandT, getStdGen)
import Control.Monad.ST (ST, runST)
import Data.Bifunctor (first)
import Data.List.Split (chunksOf)
import Data.Set (Set, fromList, member)
import Numeric.LinearAlgebra (Matrix, Normed (norm_2), maxElement, maxIndex, toRows, (><), (?))
import Numeric.LinearAlgebra.Devel (ColRange (AllCols), RowRange (Row), STMatrix, extractMatrix, freezeMatrix, readMatrix, thawMatrix, writeMatrix)

-- grid of the game
-- the current implementation relies on the grid being bounded by abyss
grid :: [[Char]]
grid =
    [ "XXXXXXXXXXXX"
    , "XS......X..X"
    , "X.......X.EX"
    , "X.......X..X"
    , "X.......X..X"
    , "X.......X..X"
    , "X..........X"
    , "XXXXXXXXXXXX"
    ]

gridWidth :: Int
gridWidth = length (head grid)

gridHeight :: Int
gridHeight = length grid

-- use a simple integer to represent the state,
-- makes it easier to index the Q matrix later on
newtype State = State {unState :: Int} deriving (Show, Eq, Ord)

toCoords :: State -> (Int, Int)
toCoords = (`divMod` gridWidth) . unState

fromCoords :: (Int, Int) -> State
fromCoords (row, col) = State $ row * gridWidth + col

indexedGrid :: [((Int, Int), Char)]
indexedGrid = concat $ zipWith (\rowIx row -> zipWith (\colIx value -> ((rowIx, colIx), value)) [0 ..] row) [0 ..] grid

startState :: State
startState = fromCoords . head $ fst <$> filter ((== 'S') . snd) indexedGrid

endState :: State
endState = fromCoords . head $ fst <$> filter ((== 'E') . snd) indexedGrid

abyss :: Set State
abyss = fromList $ fromCoords . fst <$> filter ((== 'X') . snd) indexedGrid

data Action = U | D | L | R deriving (Show, Enum)

actionToInt :: Action -> Int
actionToInt = fromEnum

intToAction :: Int -> Action
intToAction = toEnum

nextState :: State -> Action -> State
nextState state action =
    case action of
        U -> fromCoords (row - 1, col)
        D -> fromCoords (row + 1, col)
        L -> fromCoords (row, col - 1)
        R -> fromCoords (row, col + 1)
  where
    (row, col) = toCoords state

reward :: State -> Double
reward state
    | isDead state = -100
    | isGoal state = 100
    | otherwise = -1

isDead :: State -> Bool
isDead state
    | row == 0 = True
    | row == gridHeight - 1 = True
    | col == 0 = True
    | col == gridWidth - 1 = True
    | state `member` abyss = True
    | otherwise = False
  where
    (row, col) = toCoords state

isGoal :: State -> Bool
isGoal = (== endState)

-- pure type of the Q matrix
type Q = Matrix Double

-- most operations on the Q matrix will be performed
-- in the `ST` monad to avoid unnecessary copies
type STQ s = STMatrix s Double

initialQ :: Q
initialQ = (gridHeight * gridWidth) >< 4 $ repeat 0

randomAction :: (MonadRandom m) => m Action
randomAction = do
    actionAsInt <- getRandomR (0, 3)
    pure $ intToAction actionAsInt

epsilonGreedyAction :: (RandomGen g) => Double -> STQ s -> State -> RandT g (ST s) Action
epsilonGreedyAction epsilon q state = do
    randomNumber <- getRandomR (0, 1)
    if randomNumber < epsilon
        then randomAction
        else lift $ bestActionST q state

bestAction :: Q -> State -> Action
bestAction q (State s) = intToAction $ snd $ maxIndex rewards
  where
    rewards = q ? [s]

bestActionST :: STQ s -> State -> ST s Action
bestActionST q (State s) = do
    rewards <- extractMatrix q (Row s) AllCols
    pure $ intToAction $ snd $ maxIndex rewards

runGame :: Q -> State -> Double -> Int -> Maybe (State, Double)
runGame q initialState lambda maxTurns = go initialState 0 0
  where
    go state cumulativeReward turns
        | turns > maxTurns = Nothing
        | otherwise =
            let action = bestAction q state
                state' = nextState state action
                cumulativeReward' = cumulativeReward + (lambda ** fromIntegral turns) * reward state'
             in if isGoal state' || isDead state'
                    then Just (state', cumulativeReward')
                    else go state' cumulativeReward' (turns + 1)

nextQ :: Int -> Double -> State -> Action -> State -> STQ s -> ST s ()
nextQ step lambda (State s) action state'@(State s') q = do
    curVal <- readMatrix q s actionIx
    nextRow <- extractMatrix q (Row s') AllCols
    let nextVal = (1 - alpha) * curVal + alpha * r + alpha * lambda * maxElement nextRow
    writeMatrix q s actionIx nextVal
  where
    actionIx = actionToInt action
    r = reward state'
    alpha = 0.5 / fromIntegral (step + 1)

trainEpisode :: (RandomGen g) => Int -> Double -> Double -> State -> STQ s -> RandT g (ST s) ()
trainEpisode maxSteps epsilon lambda initialState q =
    do
        go 0 initialState q
  where
    go step state qST
        | step > maxSteps = pure ()
        | otherwise = do
            action <- epsilonGreedyAction epsilon qST state
            let state' = nextState state action
            lift $ nextQ step lambda state action state' qST
            if isDead state' || isGoal state'
                then pure ()
                else go (step + 1) state' qST

train :: (RandomGen g) => g -> Int -> Int -> Double -> Double -> Double -> State -> Q -> (Q, Int)
train g maxEpisodes maxSteps converged epsilon lambda state q0 =
    runST $ do
        qST <- thawMatrix q0
        episodes <- evalRandT (go 0 q0 qST) g
        q' <- freezeMatrix qST
        pure (q', episodes)
  where
    go episode qPrev qST
        | episode > maxEpisodes = pure maxEpisodes
        | otherwise =
            do
                trainEpisode maxSteps epsilon lambda state qST
                if episode > 0 && episode `mod` 100 == 0
                    then do
                        q <- lift $ freezeMatrix qST
                        if norm_2 (q - qPrev) < converged
                            then pure episode
                            else go (episode + 1) q qST
                    else go (episode + 1) qPrev qST

printQ :: Q -> String
printQ = unlines . chunksOf gridWidth . fmap (toDir . maxIndex) . toRows
  where
    toDir ix = case intToAction ix of
        U -> '^'
        D -> 'v'
        L -> '<'
        R -> '>'

main :: IO ()
main = do
    g <- getStdGen
    let q = initialQ
    let (maxEpisodes, steps, terminated, epsilon, lambda) = (1000000, 1000, 1e-9, 0.4, 0.9)
    let (trainedQ, episodes) = train g maxEpisodes steps terminated epsilon lambda startState q
    putStrLn "### training results ###\n"
    putStrLn $ "number of episodes  : " ++ show episodes
    putStrLn $ "result of tracing Q : " ++ show (first toCoords <$> runGame trainedQ startState 0.9 100)
    putStrLn $ "visualization of Q  : \n" ++ printQ trainedQ
