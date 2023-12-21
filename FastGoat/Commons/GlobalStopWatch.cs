using System.Diagnostics;

namespace FastGoat.Commons;

/// <summary>
/// A static class providing global stopwatch functionality
/// </summary>
public static class GlobalStopWatch
{
    // Stopwatch to be used by all functions in the class
    static Stopwatch sw { get; }

    // Stack to keep track of individual laps
    private static Stack<long> Lapse { get; }

    // Initializes the stopwatch and lap stack
    static GlobalStopWatch()
    {
        sw = Stopwatch.StartNew();
        Lapse = new();
    }

    /// <summary>
    /// Restarts the stopwatch and clears the lap stack
    /// </summary>
    public static void Restart()
    {
        sw.Restart();
        Lapse.Clear();
    }

    /// <summary>
    /// Records a lap on the stopwatch
    /// </summary>
    public static void AddLap()
    {
        Lapse.Push(sw.ElapsedMilliseconds);
    }

    /// <summary>
    /// Stops the stopwatch
    /// </summary>
    public static void Stop() => sw.Stop();

    /// <summary>
    /// Outputs the elapsed time since the last lap or start
    /// </summary>
    /// <param name="label">Label to be output before the time</param>
    public static void Show(string label = "")
    {
        var start = Lapse.Count == 0 ? 0 : Lapse.Pop();
        var time = TimeSpan.FromMilliseconds(sw.ElapsedMilliseconds - start);
        if (time.Hours > 0)
            Console.WriteLine($"# {label} Time:{time.Hours}h{time.Minutes}m");
        else if (time.Minutes > 0)
            Console.WriteLine($"# {label} Time:{time.Minutes}m{time.Seconds}s");
        else if (time.Seconds > 0)
            Console.WriteLine($"# {label} Time:{time.Seconds}.{time.Milliseconds:000}s");
        else
            Console.WriteLine($"# {label} Time:{sw.ElapsedMilliseconds - start}ms");
    }

    /// <summary>
    /// Runs an action and returns the elapsed time
    /// </summary>
    /// <param name="label">Label to be output before and after the action</param>
    /// <param name="action">Action to be executed and timed</param>
    /// <param name="lbl">Flag to indicate whether to output the label or not</param>
    /// <returns>The elapsed time in milliseconds</returns>
    public static long Time(string label, Action action, bool lbl = true)
    {
        if (lbl)
            Console.WriteLine($"# {label} Start");

        sw.Restart();
        action();
        sw.Stop();

        if (lbl)
            Console.WriteLine($"# {label} Time:{sw.ElapsedMilliseconds} ms");

        return sw.ElapsedMilliseconds;
    }

    /// <summary>
    /// Runs an action multiple times and outputs the average and deviation
    /// </summary>
    /// <param name="nb">Number of times to run the action</param>
    /// <param name="label">Label to be output before and after the action</param>
    /// <param name="action">Action to be executed and timed</param>
    public static void Bench(int nb, string label, Action action)
    {
        if (nb < 1)
            throw new();

        var list = new List<long>();
        for (int i = 0; i < nb; i++)
            list.Add(Time(label, action, false));

        var avg = list.Average();
        var dev = Double.Sqrt(list.Select(t => Double.Pow(t - avg, 2)).Average());
        Console.WriteLine($"# {label} Avg Time:{(long)avg} ms Dev:{dev:F}");
    }

    // Counter for infinite loop breaker
    private static int ct = 0;

    /// <summary>
    /// Resets the infinite loop breaker counter
    /// </summary>
    public static void InfiniteLoopBreakerReset() => ct = 0;

    /// <summary>
    /// Throws an exception if a specified number of iterations have been exceeded
    /// </summary>
    /// <param name="n">maximum number of iterations before the loop is considered infinite</param>
    /// <param name="msg">custom error message displayed in case the loop is broken (optional)</param>
    public static void InfiniteLoopBreaker(int n, string msg = "")
    {
        ++ct;
        if (ct > n)
            throw new ArgumentException($"################## Infinite Loop breaker ### {msg}");
    }
}