namespace FastGoat.Commons;

public class Lipsum
{
    public static string Text15Paragraphes => ByParagraphes.Glue("\n\n");
    public static string[] ByParagraphes => Lipsum15Paragraphes.Select(l => l.Glue(" ")).ToArray();
    public static string[] BySentences => Lipsum15Paragraphes.SelectMany(l => l).ToArray();

    public static string[][] Lipsum15Paragraphes = new string[][]
        {
            new string[]
            {
                "Lorem ipsum dolor sit amet, consectetur adipiscing elit.",
                "Proin consectetur fringilla turpis.",
                "Nunc facilisis enim sit amet libero pharetra finibus.",
                "Aenean quis nibh scelerisque, venenatis nulla non, viverra diam.",
                "Donec bibendum molestie iaculis.",
                "Nulla nec dolor nec neque aliquet venenatis eu at ipsum.",
                "Donec turpis erat, pretium quis aliquet quis, vestibulum vitae diam.",
                "Ut sit amet nunc ut tellus tempus ullamcorper.",
                "Mauris consequat tincidunt nibh, vitae aliquet neque.",
                "Aliquam mattis, tellus euismod hendrerit dictum, lorem sapien porta magna, eu egestas lacus nisl non lorem.",
                "Quisque vel metus congue, bibendum lorem at, feugiat nisi.",
                "Donec vestibulum dapibus dapibus.",
                "Maecenas nec nisi eu purus vehicula scelerisque.",
                "Quisque ac porttitor turpis.",
            },
            new string[]
            {
                "Duis sed ante nunc.",
                "Etiam in feugiat augue.",
                "Praesent eget auctor augue.",
                "Fusce in suscipit nunc.",
                "Mauris commodo ultricies libero sit amet eleifend.",
                "Sed quam leo, posuere at mauris in, aliquet eleifend neque.",
                "Suspendisse nec luctus orci.",
                "Sed fermentum id urna non molestie.",
                "Pellentesque vel dui molestie, ornare purus sit amet, laoreet arcu.",
                "Phasellus molestie lacinia felis at ornare.",
                "Integer lacinia sapien urna, eu lobortis nisl imperdiet a.",
                "Proin eget ante venenatis, sollicitudin mi in, tincidunt elit.",
                "Curabitur sodales euismod consectetur.",
                "Sed quis dui elit.",
                "Vivamus nisl turpis, fermentum at ultrices non, varius at mauris.",
                "Phasellus egestas ultricies dolor non fringilla.",
            },
            new string[]
            {
                "Nulla a sem luctus, dapibus sapien at, pulvinar purus.",
                "Praesent scelerisque dolor id magna blandit ultrices.",
                "Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus.",
                "Quisque non cursus ante, a consectetur urna.",
                "Sed sodales dapibus gravida.",
                "Phasellus ac dui vestibulum purus ullamcorper ornare.",
                "Etiam urna nisl, eleifend nec blandit eu, hendrerit sit amet lacus.",
                "Mauris placerat sagittis molestie.",
                "Nullam sed venenatis ante.",
            },
            new string[]
            {
                "Nulla mattis, massa sit amet aliquam pretium, velit orci malesuada arcu, nec cursus velit nisl non risus.",
                "Integer tellus sapien, sagittis nec tincidunt a, dictum nec leo.",
                "Ut quis elit id augue pretium rutrum iaculis sed justo.",
                "Cras lacus nibh, sodales ac rutrum vel, pretium convallis dui.",
                "Ut at ante sapien.",
                "Integer est nisl, hendrerit vitae mollis at, ultricies nec ipsum.",
                "In sagittis ornare dolor, vitae hendrerit lectus fringilla id.",
                "Suspendisse pharetra nisl et elementum semper.",
                "Donec fringilla, massa vitae fringilla varius, velit metus lacinia tortor, id congue mi turpis non erat.",
                "Pellentesque arcu nunc, laoreet condimentum hendrerit et, pulvinar dictum orci.",
                "Cras ante metus, pulvinar ut leo in, semper convallis sem.",
            },
            new string[]
            {
                "Ut eu sagittis tellus, ut rhoncus sapien.",
                "Duis commodo, lectus eget ultrices sagittis, eros purus ornare ligula, non varius ligula felis at eros.",
                "Fusce euismod elementum eros.",
                "Ut non nisl tortor.",
                "Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus.",
                "Sed id felis ultricies, iaculis erat eu, vehicula massa.",
                "Duis tincidunt magna nec tempus ullamcorper.",
                "Sed sed vulputate mi.",
                "Duis pharetra quam consectetur, suscipit risus sit amet, porttitor metus.",
                "Suspendisse vitae neque sed odio tincidunt tincidunt.",
                "Nam ac augue nunc.",
            },
            new string[]
            {
                "Vivamus condimentum fringilla nunc, blandit egestas purus dapibus vitae.",
                "Fusce fringilla nisi nunc, a eleifend nulla malesuada eget.",
                "Aliquam vestibulum aliquet lectus.",
                "Maecenas id ornare tellus.",
                "Ut quis sapien leo.",
                "Donec ultricies, arcu eget commodo ornare, massa est sollicitudin elit, id rhoncus diam lorem id lorem.",
                "Duis scelerisque, ex ut molestie vestibulum, dui mi condimentum elit, nec pellentesque nulla leo at nulla.",
                "Curabitur eget ultricies est, eget euismod nunc.",
                "In interdum, ante et feugiat elementum, lectus risus tempor urna, eu placerat arcu metus a enim.",
                "Cras vulputate nulla non arcu pharetra, vitae molestie nunc interdum.",
                "Vestibulum venenatis erat sed diam fringilla maximus.",
                "Donec tempus tincidunt sapien aliquet scelerisque.",
                "Proin a dignissim leo.",
            },
            new string[]
            {
                "Maecenas sagittis efficitur sem, vel venenatis ipsum.",
                "Aenean quis est a nisi placerat vulputate in pharetra felis.",
                "Sed gravida lacus at orci dignissim, eget dapibus ex molestie.",
                "Sed libero turpis, ultricies eu tortor quis, aliquet pellentesque nisi.",
                "Curabitur at pulvinar erat, eget hendrerit leo.",
                "Nulla auctor, tortor at aliquet eleifend, leo metus volutpat nulla, id lacinia leo ligula a augue.",
                "Curabitur dapibus finibus aliquam.",
                "Etiam et enim convallis, convallis odio sed, consequat lectus.",
            },
            new string[]
            {
                "Cras eget neque gravida, placerat justo fringilla, scelerisque mi.",
                "Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas.",
                "Cras mattis, leo sed finibus cursus, erat mi malesuada lorem, non congue velit quam non urna.",
                "Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Maecenas laoreet viverra lorem, a malesuada lorem vehicula ac.",
                "Donec imperdiet leo vel sem suscipit tincidunt.",
                "Cras elit dolor, faucibus vel mollis eget, condimentum ac felis.",
                "Fusce pellentesque, libero non fermentum porta, urna turpis laoreet massa, a suscipit ante nulla eget leo.",
                "Ut quis varius libero.",
                "Curabitur ipsum nisi, sollicitudin vel mi vitae, mattis ullamcorper arcu.",
                "Mauris metus magna, dignissim et mattis imperdiet, placerat sed eros.",
                "Quisque consectetur ornare aliquet.",
                "Nam aliquam ante vel pulvinar cursus.",
            },
            new string[]
            {
                "Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Proin elementum velit vel lectus ornare auctor.",
                "Cras semper magna accumsan sodales malesuada.",
                "Vestibulum ut nisi sed lorem faucibus vestibulum ut eu purus.",
                "Maecenas aliquet diam at neque venenatis dignissim.",
                "Phasellus lobortis sodales dui vel hendrerit.",
                "Nullam id nibh at arcu pulvinar tempus non sed leo.",
                "Suspendisse luctus viverra blandit.",
                "Curabitur lobortis dignissim libero nec dignissim.",
                "Ut nec tellus volutpat, fermentum metus a, iaculis metus.",
                "Curabitur eu ligula sed justo interdum pretium.",
                "Mauris nec leo facilisis, imperdiet enim at, accumsan augue.",
                "Cras non blandit ante, sit amet luctus risus.",
                "Ut orci purus, varius id placerat vitae, varius sed mauris.",
                "Aenean id mollis turpis, eu condimentum lacus.",
                "Suspendisse auctor arcu a eros tristique, et vestibulum ligula varius.",
            },
            new string[]
            {
                "Vivamus vel orci finibus, aliquam neque quis, porttitor quam.",
                "Proin scelerisque, lorem id tincidunt malesuada, lacus metus eleifend mi, ut condimentum est tellus ut nisi.",
                "Nulla facilisi.",
                "In molestie nunc id ultricies consequat.",
                "Pellentesque tristique velit vulputate massa mollis, non laoreet orci venenatis.",
                "Aenean in diam vitae purus congue luctus.",
                "Sed hendrerit rutrum tortor non convallis.",
                "In pharetra lorem enim, vitae euismod turpis viverra sed.",
                "Nunc nec ligula ut quam iaculis venenatis quis quis massa.",
                "In fermentum risus nec sem porta, in dictum erat sollicitudin.",
            },
            new string[]
            {
                "Nam viverra est vel quam euismod, in ultrices est ultrices.",
                "Duis sed leo augue.",
                "Morbi id felis vestibulum, finibus tortor sed, molestie lacus.",
                "Ut malesuada dolor ut ex congue rutrum.",
                "Ut iaculis libero ante, sed dictum nibh ultricies quis.",
                "Sed mattis et sem at volutpat.",
                "Pellentesque ultricies volutpat pharetra.",
                "Proin sit amet magna est.",
            },
            new string[]
            {
                "Duis ac mi ut sapien rutrum accumsan vitae nec felis.",
                "Duis purus quam, condimentum id tellus commodo, semper ultrices urna.",
                "Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus.",
                "Nulla dictum nisi massa, vitae porttitor ipsum pulvinar a.",
                "Donec vitae nisi maximus, aliquet sem sed, consectetur est.",
                "Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos.",
                "Duis posuere risus nec ultrices facilisis.",
                "Suspendisse aliquam urna ac sodales viverra.",
                "Curabitur tincidunt erat in quam dignissim, quis elementum purus ornare.",
                "Suspendisse ut porttitor nisi.",
                "Suspendisse interdum lectus malesuada ligula interdum, in iaculis ex pharetra.",
                "Maecenas luctus pharetra gravida.",
            },
            new string[]
            {
                "Duis ac felis eget nibh scelerisque vulputate at id dui.",
                "Proin finibus non nisl non efficitur.",
                "Nulla at lacus in mi tristique pretium.",
                "Mauris hendrerit erat turpis, vitae vehicula tortor euismod ac.",
                "Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas.",
                "Nulla facilisi.",
                "Vestibulum pharetra ex ante, aliquam consectetur felis blandit nec.",
                "Suspendisse non quam in justo cursus facilisis eget sit amet urna.",
                "Sed id mattis nunc.",
                "Proin in elit nec quam euismod molestie.",
                "Maecenas euismod, libero ac finibus viverra, magna justo iaculis ligula, ut commodo magna risus vitae eros.",
                "Morbi fermentum ante felis, vel sodales dui pharetra in.",
                "Suspendisse volutpat laoreet varius.",
                "Aliquam condimentum eros ut laoreet mollis.",
                "Suspendisse potenti.",
                "Nullam vulputate enim diam, finibus lobortis sapien mollis nec.",
            },
            new string[]
            {
                "Phasellus tellus tellus, molestie id blandit eget, consequat eget tortor.",
                "Mauris dui dolor, lacinia quis justo ut, consectetur porttitor odio.",
                "Phasellus eget fringilla turpis, id venenatis magna.",
                "Quisque porttitor velit massa, eget lobortis purus mollis nec.",
                "Morbi cursus pretium posuere.",
                "Nunc nec nulla orci.",
                "Pellentesque eget consequat sapien.",
                "Nunc semper vulputate libero sit amet blandit.",
            },
            new string[]
            {
                "Donec eget arcu a velit egestas volutpat.",
                "Proin justo justo, tristique ac convallis sit amet, suscipit a metus.",
                "Donec nec auctor tellus.",
                "Nulla commodo nisi condimentum erat vestibulum, quis suscipit libero bibendum.",
                "Quisque varius sodales sodales.",
                "Praesent massa purus, egestas quis ante ut, ullamcorper vehicula augue.",
                "Fusce eget bibendum nisl.",
            },
        };
}